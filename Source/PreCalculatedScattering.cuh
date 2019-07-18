/*
Copyright (c) 2081, W. Groen <w.a.groun@student.tudelft.nl>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include "Transport.cuh"
#include <curand_kernel.h>

#define OPACITY_BLOCK_DIM	8

KERNEL void KrnlPreCalculatedScattering(CScene* pScene, CCudaView* pView, Vec3f* pPoints, int nrPoints, float* pDevConnections, CColorXyza* pDevPointColour, Vec3f* pResults)
{
	// Get X and Y position in the block grid.
	const int X = blockIdx.x * blockDim.x + threadIdx.x;
	const int Y = blockIdx.y * blockDim.y + threadIdx.y;

	// Check if within bounds of the canvas
	if (X >= gFilmWidth || Y >= gFilmHeight)
		return;

	// Create RNG using 2 random seeds
	CRNG RNG(pView->m_RandomSeeds1.GetPtr(X, Y), pView->m_RandomSeeds2.GetPtr(X, Y));

	// Define base colours
	// Lv = final colour
	// Li = colour of currently selected light, gets filled during certain calls
	//CColorXyz Lv = SPEC_BLACK,
	CColorXyz Li = SPEC_CYAN;// SPEC_BLACK;
	CColorXyza Lv = CColorXyza(0, 1, 0);

	// Declare the ray
	CRay Re;

	// Define the pixel
	// UV = current pixel in grid + ?random offset?
	const Vec2f UV = Vec2f(X, Y) + RNG.Get2();

	// Set the ray's origin and destination
	// re.m_O is the position of the camera
	// re.m_D is the direction the camera is facing
	// If the aperture of the camere is greater than 0, the RNG factor helps add some blur to the background
	pScene->m_Camera.GenerateRay(UV, RNG.Get2(), Re.m_O, Re.m_D);

	//CRay Re2;
	//CRay Re3;
	//CRay Re4;
	//const Vec2f UV2 = Vec2f(X, Y + 10);
	//const Vec2f UV3 = Vec2f(X + 10, Y);
	//const Vec2f UV4 = Vec2f(X + 10, Y+10);
	//pScene->m_Camera.GenerateRay(UV2, RNG.Get2(), Re2.m_O, Re2.m_D);
	//pScene->m_Camera.GenerateRay(UV3, RNG.Get2(), Re3.m_O, Re3.m_D);
	//pScene->m_Camera.GenerateRay(UV4, RNG.Get2(), Re4.m_O, Re4.m_D);

	//pResults[Y * gFilmWidth + X] = Vec3f((Re2.m_D - Re.m_D).Length(), (Re3.m_D - Re.m_D).Length(), (Re4.m_D - Re.m_D).Length());
	//pResults[Y * gFilmWidth + X] = Re.m_D;
	//return;

	// Set Woodcock info?????????
	// m_MinT = ??
	// m_MaxT = ??
	Re.m_MinT = 0.0f;
	Re.m_MaxT = 1500000.0f;

	// Declare ???????
	// Pe = The position of the scattering point. Gets filled by SampleDistance()
	// Pl = ?????
	Vec3f Pe, Pl;

	// Declare ???????
	CLight* pLight = NULL;

	// Check if we hit the boundingbox
	float MinT = 0;
	float MaxT = 0;
	if (!IntersectBox(Re, &MinT, &MaxT)) {
		if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, INF_MAX), Li, Pl, pLight))
			Lv = CColorXyza(Li.c[0], Li.c[1], Li.c[2]);
		else
			Lv = CColorXyza(0, 0, 1);
	}
	else {

		int closest = -1;
		float dist = 999999999; // some large number
		float threshold = 0.001;
		float minD = 99999999999999;
		for (int i = 0; i < nrPoints; i++) {
			//if (X == 32 && Y == 32)
			//	pResults[i] = pPoints[i];
			if (pDevPointColour[i].IsBlack())
				continue;

			Vec3f point = pPoints[i];
			Vec3f aMinP = point - Re.m_O;
			float length = aMinP.Dot(Re.m_D);
			Vec3f pointOnLine = (length * Re.m_D) + Re.m_O;
			float d = (point - pointOnLine).Length();

			//pResults[Y * gFilmWidth + X] = Vec3f(d, -1,-1);

			minD = Fminf(minD, d);

			if (d < threshold) {
				float distFromCamera = length;
				if (distFromCamera > 0 && distFromCamera < dist) {
					closest = i;
					dist = distFromCamera;
					//pResults[Y * gFilmWidth + X] = Vec3f(i, d, length);
					//pResults[Y * gFilmWidth + X] = point;
				}
			}

		}

		if (closest >= 0 && closest < nrPoints) {
			Lv = pDevPointColour[closest];
			//Lv = CColorXyza(0, 0, 1);
		}

		if (closest == -1) {
			//Lv = CColorXyza(minD, minD, minD);

			if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, INF_MAX), Li, Pl, pLight))
				Lv = CColorXyza(Li.c[0], Li.c[1], Li.c[2]);
			else
				Lv = CColorXyza(1, 1, 0);
		}
	}

	__syncthreads();

	pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
	//pView->m_FrameEstimateXyza.Set(CColorXyza(RNG.Get1(), RNG.Get1(), RNG.Get1()), X, Y);

	//pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
	//int t = Y * gFilmWidth + X;
	//if (t < nrPoints)
	//	pView->m_FrameEstimateXyza.Set(pDevPointColour[t], X, Y);
	//else
	//	pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
}

void PreCalculatedScattering(CScene* pScene, CScene* pDevScene, CCudaView* pView, Vec3f* pPoints, int nrPoints, float* pDevConnections, CColorXyza* pDevPointColour, Vec3f* pDevResults)
{
	const dim3 KernelBlock(KRNL_SS_BLOCK_W, KRNL_SS_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResX() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResY() / (float)KernelBlock.y));

	KrnlPreCalculatedScattering<<<KernelGrid, KernelBlock>>>(pDevScene, pView, pPoints, nrPoints, pDevConnections, pDevPointColour, pDevResults);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Pre-calculated Scattering");
}

//curandState initStates(curandState* states) {
//	return
//}

KERNEL void krnlGeneratePoints(CScene* scene, curandState* states, Vec3f* points, int nrStates, int pointsPerState, int nrPoints) {
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id >= nrStates)
		return;
	
	curandState localState = states[id];

	Vec3f base = Vec3f(scene->m_BoundingBox.LengthX(), scene->m_BoundingBox.LengthY(), scene->m_BoundingBox.LengthZ());
	for (int i = 0; i < pointsPerState; i++) {
		//float r = curand_uniform(&localState);
		int idPoint = id * pointsPerState + i;
		
		if (idPoint < nrPoints) {
			//points[idPoint] = Vec3f(base) *r;
			points[idPoint] = Vec3f(base.x * curand_uniform(&localState), base.y * curand_uniform(&localState), base.z * curand_uniform(&localState));
		}
	}

	states[id] = localState;
}

void GeneratePoints(CScene* devScene, curandState* states, Vec3f* points, int nrStates, int pointsPerState, int nrPoints) {
	int threadsPerBlock = 512;
	int blocks = (nrStates + (threadsPerBlock - 1)) / threadsPerBlock;

	krnlGeneratePoints<<<blocks, threadsPerBlock >>>(devScene, states, points, nrStates, pointsPerState, nrPoints);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Generating Points");
}

KERNEL void krnlSetupStates(curandState* states, int* seeds, int N) {
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id < N)
		curand_init(seeds[id], id, 0, &states[id]);
}

KERNEL void krnlSetupStatesSingleSeed(curandState* states, int seed, int N) {
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id < N)
		curand_init(seed, id, 0, &states[id]);
}

void SetupStates(curandState* states, int* seeds, int n) {
	int threadsPerBlock = 256;
	int blocks = (n + (threadsPerBlock - 1)) / threadsPerBlock;

	krnlSetupStates<<<blocks, threadsPerBlock >>>(states, seeds, n);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Setting up curandStates");
}

void SetupStatesSingleSeed(curandState* states, int seed, int n) {
	int threadsPerBlock = 512;
	int blocks = (n + (threadsPerBlock - 1)) / threadsPerBlock;

	krnlSetupStatesSingleSeed<<<blocks, threadsPerBlock >>>(states, seed, n);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Setting up curandStates");
}


KERNEL void krnlMakeConnections(CScene* pScene, Vec3f* pPoints, int nrPoints, unsigned int* pSeed0, unsigned int* pSeed1, float* pDevConnections, CColorXyza* pDevPointColour) {
	// Get X and Y position in the block grid.
	const int X = blockIdx.x * blockDim.x + threadIdx.x;
	const int Y = blockIdx.y * blockDim.y + threadIdx.y;

	const int id = X;

	// Start by calculating the colour of the point
	if (id < nrPoints && Y == nrPoints - 1)
	{
		CRNG RNG(&pSeed0[id], &pSeed1[id]);

		Vec3f Pe, Pl;
		Pe = pPoints[id];


		CColorXyz Lv = SPEC_BLACK, Li = SPEC_BLACK;
		CLight* pLight = NULL;

		const float D = GetNormalizedIntensity(Pe);
		Lv += GetEmission(D).ToXYZ();

		// I wonder how much influence the direction of the ray holds
		Vec3f Re = Vec3f(0.5, 0.5, 0.5);

		// Switch Depending on the shading type
		switch (pScene->m_ShadingType)
		{
			// BRDF Only (Bidirectional Reflectance Distribution Function)
			case 0:
			{
				Lv += UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re), Pe, NormalizedGradient(Pe), RNG, true);
				break;
			}

			// Phase Function Only
			case 1:
			{
				Lv += 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re), Pe, NormalizedGradient(Pe), RNG, false);
				break;
			}

			// Hybrid (BDRF & Phase Function)
			case 2:
			{

				const float GradMag = GradientMagnitude(Pe) * gIntensityInvRange;

				const float PdfBrdf = (1.0f - __expf(-pScene->m_GradientFactor * GradMag));

				if (RNG.Get1() < PdfBrdf)
					Lv += UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re), Pe, NormalizedGradient(Pe), RNG, true);
				else
					Lv += 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re), Pe, NormalizedGradient(Pe), RNG, false);

				break;
			}
		}

		pDevPointColour[id] = CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]);
	}
	return;
	// calculate throughput to other points
	if (X < nrPoints && Y < nrPoints && Y > X)
	{
		//for (int i = id + 1; i < nrPoints; i++) {
		Vec3f Pe = pPoints[X];
		Vec3f Pe2 = pPoints[Y];

		// TODO get actual value
		float blockage = 1 - (1 / (Pe2 - Pe).Length());

		//int id2 = id * nrPoints + i;
		pDevConnections[X] = blockage;
		//int id2_2 = i * nrPoints + id;
		pDevConnections[Y] = blockage;

		//pDevConnections[id].resize(nrPoints);
		//pDevConnections[id].at(i) = blockage;

		//pDevConnections[i].resize(nrPoints);
		//pDevConnections[i].at(id) = blockage;
		//}
	}
}

void makeConnections(CScene* pDevScene, Vec3f* pDevPoints, int nrPoints, unsigned int* pDevSeed0, unsigned int* pDevSeed1, float* pDevConnections, CColorXyza* pDevPointColour)
{
	//int threadsPerBlock = 256;
	//int blocks = (nrPoints + (threadsPerBlock - 1)) / threadsPerBlock;


	const dim3 KernelBlock(KRNL_SS_BLOCK_W, KRNL_SS_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)nrPoints / (float)KernelBlock.x), (int)ceilf((float)nrPoints / (float)KernelBlock.y));

	//krnlMakeConnections<<<blocks, threadsPerBlock>>>(pDevScene, pDevPoints, nrPoints, pDevSeed0, pDevSeed1, pDevConnections, pDevPointColour);
	krnlMakeConnections<<<KernelGrid, KernelBlock >>>(pDevScene, pDevPoints, nrPoints, pDevSeed0, pDevSeed1, pDevConnections, pDevPointColour);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Calculating connections");
}

KERNEL void KrnlGetScatterInfo(CScene* pScene, CCudaView* pView, Vec3f* pPoints, curandState* pStates)
{
	// Get X and Y position in the block grid.
	const int X = blockIdx.x * blockDim.x + threadIdx.x;
	const int Y = blockIdx.y * blockDim.y + threadIdx.y;

	// Check if within bounds of the canvas
	if (X >= gFilmWidth || Y >= gFilmHeight)
		return;

	// Create RNG using 2 random seeds
	CRNG RNG(pView->m_RandomSeeds1.GetPtr(X, Y), pView->m_RandomSeeds2.GetPtr(X, Y));

	int stateId = Y * gFilmWidth + X;
	curandState localState = pStates[stateId];

	//float r = curand_uniform(&localState);
	//float r = RNG.Get1();
	//pView->m_FrameEstimateXyza.Set(CColorXyza(r, r, r), X, Y);
	//pPoints[stateId] = Vec3f(r, r, r);

	//pStates[stateId] = localState;

	//return;



	CColorXyz Lv = SPEC_BLACK, Li = SPEC_BLACK;

	// Declare the ray
	CRay Re;

	// Define the pixel
	// UV = current pixel in grid + ?random offset?
	const Vec2f UV = (pScene->m_PostProcessingSteps & 16) ? Vec2f(X, Y) + RNG.Get2() : Vec2f(X, Y);
	//const Vec2f UV = (pScene->m_PostProcessingSteps & 16) ? Vec2f(X, Y) + Vec2f(curand_uniform(&state), curand_uniform(&state)) : Vec2f(X, Y);

	// Set the ray's origin and destination
	// re.m_O is the position of the camera
	// re.m_D is the direction the camera is facing
	// If the aperture of the camere is greater than 0, the RNG factor helps add some blur to the background
	pScene->m_Camera.GenerateRay(UV, RNG.Get2(), Re.m_O, Re.m_D);

	// Set Woodcock info?????????
	// m_MinT = ??
	// m_MaxT = ??
	Re.m_MinT = 0.0f;
	Re.m_MaxT = 1500000.0f;

	// Declare ???????
	// Pe = The position of the scattering point. Gets filled by SampleDistance()
	// Pl = ?????
	Vec3f Pe, Pl;

	// Declare ???????
	CLight* pLight = NULL;

	// Find a suitable scattering point. If not return false. Position is stored in Pe
	if (SampleDistanceRM(Re, RNG, Pe))
	{
		pPoints[Y * gFilmWidth + X] = Pe;
	}
	//pPoints[Y * gFilmWidth + X] = RNG.Get3();
}

void getScatterInfo(CScene* pScene, CScene* pDevScene, CCudaView* pView, Vec3f* pPoints, curandState* pStates) {
	const dim3 KernelBlock(KRNL_SS_BLOCK_W, KRNL_SS_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResX() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResY() / (float)KernelBlock.y));
	// Some weird block distortion in a 17x21 grid at 400x500 pixels, scales with the amount of pixels
	KrnlGetScatterInfo <<<KernelGrid, KernelBlock>>>(pDevScene, pView, pPoints, pStates);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Single Scattering Point Info");
}

KERNEL void KrnlCreateOpacityGradientTexture(CResolution3D Resolution, float4* pDevOpacityGradient1D, float* pDevOpacityGradientMagnitude1D) {
	const int X = blockIdx.x * blockDim.x + threadIdx.x;
	const int Y = blockIdx.y * blockDim.y + threadIdx.y;
	const int Z = blockIdx.z * blockDim.z + threadIdx.z;

	if (X >= Resolution.GetResX() || Y >= Resolution.GetResY() || Z >= Resolution.GetResZ())
		return;

	__shared__ float opacityStored[OPACITY_BLOCK_DIM + 2][OPACITY_BLOCK_DIM + 2][OPACITY_BLOCK_DIM + 2];

	Vec3f realPos = gGradientDelta * Vec3f(X, Y, Z);

	// Position in shared memory
	int memoryX = threadIdx.x + 1;
	int memoryY = threadIdx.y + 1;
	int memoryZ = threadIdx.z + 1;

	opacityStored[memoryX][memoryY][memoryZ] = GetOpacity(GetNormalizedIntensity(realPos));

	// Get the extra points for the x edges
	if (memoryX == 1) {
		opacityStored[0][memoryY][memoryZ] = GetOpacity(GetNormalizedIntensity(realPos - ToVec3f(gGradientDeltaX)));
	}
	else if (memoryX == blockDim.x) {
		opacityStored[memoryX + 1][memoryY][memoryZ] = GetOpacity(GetNormalizedIntensity(realPos + ToVec3f(gGradientDeltaX)));
	}

	// Get the extra points for the y edges
	if (memoryY == 1) {
		opacityStored[memoryX][0][memoryZ] = GetOpacity(GetNormalizedIntensity(realPos - ToVec3f(gGradientDeltaY)));
	}
	else if (memoryY == blockDim.y) {
		opacityStored[memoryX][memoryY + 1][memoryZ] = GetOpacity(GetNormalizedIntensity(realPos + ToVec3f(gGradientDeltaY)));
	}

	// Get the extra points for the z edges
	if (memoryZ == 1) {
		opacityStored[memoryX][memoryY][0] = GetOpacity(GetNormalizedIntensity(realPos - ToVec3f(gGradientDeltaZ)));
	}
	else if (memoryZ == blockDim.z) {
		opacityStored[memoryX][memoryY][memoryZ + 1] = GetOpacity(GetNormalizedIntensity(realPos + ToVec3f(gGradientDeltaZ)));
	}

	// Sync all threads to ensure memory is filled
	__syncthreads();

	// Calculate the gradient
	float4 gradient;
	gradient.x = opacityStored[memoryX + 1][memoryY][memoryZ] - opacityStored[memoryX - 1][memoryY][memoryZ];
	gradient.y = opacityStored[memoryX][memoryY + 1][memoryZ] - opacityStored[memoryX][memoryY - 1][memoryZ];
	gradient.z = opacityStored[memoryX][memoryY][memoryZ + 1] - opacityStored[memoryX][memoryY][memoryZ - 1];

	// Write the gradient to the 1D array
	pDevOpacityGradient1D[Index3To1(X, Y, Z, Resolution)] = gradient;
	pDevOpacityGradientMagnitude1D[Index3To1(X, Y, Z, Resolution)] = Vec3f(gradient.x, gradient.y, gradient.z).Length();
}

void GetOpacityGradientTexture(CResolution3D resolution, float4* pDevOpacityGradient1D, float* pDevOpacityGradientMagnitude1D) {
	//const dim3 KernelBlock(512, 512, 512);
	//const dim3 KernelGrid((int)ceilf((float)Resolution.GetResX() / (float)KernelBlock.x), (int)ceilf((float)Resolution.GetResY() / (float)KernelBlock.y), (int)ceilf((float)Resolution.GetResZ() / (float)KernelBlock.z));

	const dim3 KernelBlock((int)ceilf((float)resolution.GetResX() / (float)OPACITY_BLOCK_DIM), (int)ceilf((float)resolution.GetResY() / (float)OPACITY_BLOCK_DIM), (int)ceilf((float)resolution.GetResZ() / (float)OPACITY_BLOCK_DIM));
	const dim3 KernelGrid(OPACITY_BLOCK_DIM, OPACITY_BLOCK_DIM, OPACITY_BLOCK_DIM);
	
	KrnlCreateOpacityGradientTexture<<<KernelBlock, KernelGrid>>>(resolution, pDevOpacityGradient1D, pDevOpacityGradientMagnitude1D);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Opacity Gradient Creation");
}


//KERNEL void KrnlCreateIlluminationTexture(CScene* pScene, float* pDevIllumination1D, CCudaRandomBuffer2D* rng1, CCudaRandomBuffer2D* rng2) {
KERNEL void KrnlCreateIlluminationTexture(CScene* pScene, float* pDevIllumination1D, unsigned int* rng1, unsigned int* rng2) {
	const int X = blockIdx.x * blockDim.x + threadIdx.x;
	const int Y = blockIdx.y * blockDim.y + threadIdx.y;
	const int Z = blockIdx.z * blockDim.z + threadIdx.z;

	if (X >= pScene->m_Resolution.GetResX() || Y >= pScene->m_Resolution.GetResY() || Z >= pScene->m_Resolution.GetResZ())
		return;

	Vec3f realPos = gGradientDelta * Vec3f(X, Y, Z);

	//__shared__ CRNG rngTest[OPACITY_BLOCK_DIM][OPACITY_BLOCK_DIM];
	//rngTest[threadIdx.x][threadIdx.y] = CRNG(rng1->GetPtr(X, Y), rng2->GetPtr(X, Y));
	//CRNG RNG = rngTest[threadIdx.x][threadIdx.y];
	//CRNG RNG = CRNG(rng1->GetPtr(X, Y), rng2->GetPtr(X, Y));
	CRNG RNG = CRNG(&rng1[Index3To1(X, Y, Z, pScene->m_Resolution)], &rng2[Index3To1(X, Y, Z, pScene->m_Resolution)]);
	
	CColorXyz Lv = SPEC_BLACK, Li = SPEC_BLACK, Tr = SPEC_WHITE;
	
	CRay Re;

	Re.m_O = realPos;
	Re.m_D = UniformSampleSphere(RNG.Get2());
	Re.m_MinT = 0.0f;
	Re.m_MaxT = FLT_MAX;

	Vec3f Pe, Pl;

	CLight* pLight = NULL;
	for (int i = 0; i < pScene->m_MaxBounces; i++) {
		if (SampleDistanceRM(Re, RNG, Pe)) {
			const float D = GetNormalizedIntensity(Pe);

			Lv += Tr * GetEmission(D).ToXYZ();

			// Switch Depending on the shading type
			switch (pScene->m_ShadingType) {
				// BRDF Only (Bidirectional Reflectance Distribution Function)
				case 0:	{
					Lv += Tr * UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
					break;
				}

				// Phase Function Only
				case 1:	{
					Lv += Tr * 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);
					break;
				}

				// Hybrid (BDRF & Phase Function)
				case 2:	{
					const float GradMag = GradientMagnitude(Pe) * gIntensityInvRange;

					const float PdfBrdf = (1.0f - __expf(-pScene->m_GradientFactor * GradMag));

					if (RNG.Get1() < PdfBrdf)
						Lv += Tr * UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
					else
						Lv += Tr * 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);

					break;
				}
			}
		}
		else {
			if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, INF_MAX), Li, Pl, pLight)) {
				Lv += Tr * Li;
			}
			break;
		}

		Re.m_O = Pe;
		Re.m_D = UniformSampleSphere(RNG.Get2());
		Re.m_MinT = 0.0f;
		Re.m_MaxT = INF_MAX;

		// Adjusting weight for sampling from a sphere
		Tr *= INV_4_PI_F;

		if (Terminate(Tr, RNG)) {
			break;
		}

		
	}
	
	// Convert to grayscale
	pDevIllumination1D[Index3To1(X, Y, Z, pScene->m_Resolution)] = (Lv.c[0] + Lv.c[1] + Lv.c[2]) / 3.0f;
}

//void GetIlluminationTexture(CScene* pScene, CResolution3D resolution, float* pDevIllumination1D, CCudaRandomBuffer2D* rng1, CCudaRandomBuffer2D* rng2) {
void GetIlluminationTexture(CScene* pScene, CResolution3D resolution, float* pDevIllumination1D, unsigned int* rng1, unsigned int* rng2) {
	const dim3 KernelBlock((int)ceilf((float)resolution.GetResX() / (float)OPACITY_BLOCK_DIM), (int)ceilf((float)resolution.GetResY() / (float)OPACITY_BLOCK_DIM), (int)ceilf((float)resolution.GetResZ() / (float)OPACITY_BLOCK_DIM));
	const dim3 KernelGrid(OPACITY_BLOCK_DIM, OPACITY_BLOCK_DIM, OPACITY_BLOCK_DIM);

	KrnlCreateIlluminationTexture<<<KernelBlock, KernelGrid>>>(pScene, pDevIllumination1D, rng1, rng2);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Illumination Texture Creation");
}


KERNEL void krnlTestMultipleScattering(CScene* pScene, CCudaView* pView)
{
	const int X = (blockIdx.x * blockDim.x) + threadIdx.x;
	const int Y = (blockIdx.y * blockDim.y) + threadIdx.y;
	const int PID = (Y * gFilmWidth) + X;

	if (X >= gFilmWidth || Y >= gFilmHeight || PID >= gFilmNoPixels)
		return;

	//CRNG RNG(&pSeeds[2 * PID], &pSeeds[2 * PID + 1]);
	CRNG RNG(pView->m_RandomSeeds1.GetPtr(X, Y), pView->m_RandomSeeds2.GetPtr(X, Y));

	CColorXyz Lv = SPEC_BLACK, Li = SPEC_BLACK, Tr = SPEC_WHITE;

	CRay Re;

	const Vec2f UV = Vec2f(X, Y) + RNG.Get2();

	pScene->m_Camera.GenerateRay(UV, RNG.Get2(), Re.m_O, Re.m_D);

	Re.m_MinT = 0.0f;
	Re.m_MaxT = FLT_MAX;

	Vec3f Pe, Pl;

	CLight* pLight = NULL;
	float test = 100;
	float range = 1000;
	bool testing = false;
	for (int i = 0; i < pScene->m_MaxBounces; i++)
		//for (int i = 0; !Terminate(Tr, RNG); i++)
		//for (int i = 0; i < 1; i++)
	{
		if (SampleDistanceRM(Re, RNG, Pe))
		{
			/*if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, (Pe - Re.m_O).Length()), Li, Pl, pLight))
			{
			pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
			//float4 ColorXYZA = make_float4(Lv.c[0], Lv.c[1], Lv.c[2], 0.0f);
			//				surf2Dwrite(ColorXYZA, gSurfRunningEstimateXyza, X * sizeof(float4), Y);

			return;
			}*/

			const float D = GetNormalizedIntensity(Pe);

			Lv += Tr * GetEmission(D).ToXYZ();

			// original
			//Lv += Tr * 0.5f * UniformSampleOneLight(pScene, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);
			// fixed
			//Lv += Tr * 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);

			// from single
			// Switch Depending on the shading type
			switch (pScene->m_ShadingType)
			{
				// BRDF Only (Bidirectional Reflectance Distribution Function)
			case 0:
			{
				Lv += Tr * UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
				if (testing && (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test)) {
					pView->m_FrameEstimateXyza.Set(CColorXyza((i + 1) / (float)pScene->m_MaxBounces, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / range, 0.4f), X, Y);
					return;
				}
				break;
			}

			// Phase Function Only
			case 1:
			{
				Lv += Tr * 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);
				if (testing && (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test)) {
					pView->m_FrameEstimateXyza.Set(CColorXyza((i + 1) / (float)pScene->m_MaxBounces, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / range, 0.6f), X, Y);
					return;
				}
				break;
			}

			// Hybrid (BDRF & Phase Function)
			case 2:
			{
				const float GradMag = GradientMagnitude(Pe) * gIntensityInvRange;

				const float PdfBrdf = (1.0f - __expf(-pScene->m_GradientFactor * GradMag));
				float v = RNG.Get1() < PdfBrdf ? 0.4f : 0.6f;
				if (RNG.Get1() < PdfBrdf)
					Lv += Tr * UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
				else
					Lv += Tr * 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);

				if (testing && (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test)) {
					pView->m_FrameEstimateXyza.Set(CColorXyza((i + 1) / (float)pScene->m_MaxBounces, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / range, v), X, Y);
					return;
				}

				break;
			}
			}

			//Tr *= 1.0f - GetOpacity(D);
		}
		else
		{
			if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, INF_MAX), Li, Pl, pLight)) {
				Lv += Tr * Li;

				if (testing && (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test)) {
					pView->m_FrameEstimateXyza.Set(CColorXyza((i + 1) / (float)pScene->m_MaxBounces, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / range, 1), X, Y);
					return;
				}
			}
			if (testing && (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test)) {
				pView->m_FrameEstimateXyza.Set(CColorXyza((i + 1) / (float)pScene->m_MaxBounces, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / range, 0.8f), X, Y);
				return;
			}
			break;
		}



		Re.m_O = Pe;
		Re.m_D = UniformSampleSphere(RNG.Get2());
		Re.m_MinT = 0.0f;
		Re.m_MaxT = INF_MAX;

		// Adjusting weight for sampling from a sphere
		Tr *= INV_4_PI_F;

		if (Terminate(Tr, RNG)) {
			if (testing && (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test)) {
				pView->m_FrameEstimateXyza.Set(CColorXyza((i + 1) / (float)pScene->m_MaxBounces, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / range, 0.2f), X, Y);
				return;
			}
			break;
		}
	}

	__syncthreads();

	if (testing) {
		if (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test) {
			//pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0] / 10.0f, Lv.c[1] / 10.0f, Lv.c[2] / 10.0f), X, Y);
			//pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
			pView->m_FrameEstimateXyza.Set(CColorXyza(1, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / range, 0), X, Y);
		}
		else
			pView->m_FrameEstimateXyza.Set(CColorXyza(0), X, Y);
	}
	else {
		/*if (Lv.c[0] > test || Lv.c[1] > test || Lv.c[2] > test) {
		pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0] / 10.0f, Lv.c[1] / 10.0f, Lv.c[2] / (10.0f * test)), X, Y);
		//pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
		//pView->m_FrameEstimateXyza.Set(CColorXyza(1, max(max(Lv.c[0], Lv.c[1]), Lv.c[2]) / 10.0f, 0.25f), X, Y);
		}
		else
		pView->m_FrameEstimateXyza.Set(CColorXyza(0), X, Y);*/

		//pView->m_FrameEstimateXyza.Set(CColorXyza(Tr.c[0], Tr.c[0], Tr.c[0]), X, Y);
		pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);

		//pView->m_FrameEstimateXyza.Set(CColorXyza(Clamp(Lv.c[0], 0.0, 1.0f), Clamp(Lv.c[1], 0.0, 1.0f), Clamp(Lv.c[2], 0.0, 1.0f)), X, Y);

		//float4 ColorXYZA = make_float4(Lv.c[0], Lv.c[1], Lv.c[2], 0.0f);
		//	surf2Dwrite(ColorXYZA, gSurfRunningEstimateXyza, X * sizeof(float4), Y);
	}
}

//void MultipleScattering(CScene* pScene, CScene* pDevScene, int* pSeeds)
void TestMultipleScattering(CScene* pScene, CScene* pDevScene, CCudaView* pView)
{
	//const dim3 KernelBlock(KRNL_MS_BLOCK_W, KRNL_MS_BLOCK_H);
	const dim3 KernelBlock(KRNL_SS_BLOCK_W, KRNL_SS_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResX() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResY() / (float)KernelBlock.y));

	krnlTestMultipleScattering<<<KernelGrid, KernelBlock>>>(pDevScene, pView);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Test Multiple Scattering");
}