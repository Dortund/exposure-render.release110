/*
	Copyright (c) 2011, T. Kroes <t.kroes@tudelft.nl>
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

KERNEL void KrnlSingleScattering(CScene* pScene, CCudaView* pView)
{
	// Get X and Y position in the block grid.
	const int X		= blockIdx.x * blockDim.x + threadIdx.x;
	const int Y		= blockIdx.y * blockDim.y + threadIdx.y;

	// Check if within bounds of the canvas
	if (X >= gFilmWidth || Y >= gFilmHeight)
		return;
	
	// Create RNG using 2 random seeds
	CRNG RNG(pView->m_RandomSeeds1.GetPtr(X, Y), pView->m_RandomSeeds2.GetPtr(X, Y));

	//float r = RNG.Get1();
	//pView->m_FrameEstimateXyza.Set(CColorXyza(r, r, r), X, Y);
	//return;

	//pView->m_FrameEstimateXyza.Set(CColorXyza((float) X / (float)gFilmWidth, (float) Y / (float)gFilmHeight, fmodf(X*0.1, 1.f) ), X, Y);
	//return;

	//curandState state;
	//curand_init((int) pView->m_RandomSeeds1.GetPtr(X, Y), 1, 0, &state);
	//float r = curand_uniform(&state);

	// Define base colours
	// Lv = final colour
	// Li = colour of currently selected light, gets filled during certain calls
	CColorXyz Lv = SPEC_BLACK, Li = SPEC_BLACK;

	// Declare the ray
	CRay Re;
	
	// Define the pixel
	// UV = current pixel in grid + ?random offset?
	const Vec2f UV = (pScene->m_PostProcessingSteps & 16) ? Vec2f(X, Y) + RNG.Get2() : Vec2f(X, Y);
	//const Vec2f UV = (pScene->m_PostProcessingSteps & 16) ? Vec2f(X, Y) + Vec2f(curand_uniform(&localState), curand_uniform(&localState)) : Vec2f(X, Y);

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
	//if (SampleDistanceRM(Re, localState, Pe))
	if (SampleDistanceRM(Re, RNG, Pe))
	{
		// Retrieves information about a light it had a hit with
		// Li = light color, gets filled during the call
		// Pl = the position of the light along the ray, gets filled during the call
		// pLight = the light with which we have a hit, gets filled during the call
		if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, (Pe - Re.m_O).Length()), Li, Pl, pLight))
		{
			// Set the current estimated colour of the pixel. Seems to always be 'SPEC_BLACK' from this position in the code.
			// And return on the function call
			pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
			return;
		}
		
		// Get the normalized intensity at the scatteringpoint
		const float D = GetNormalizedIntensity(Pe);

		// Get the colour emission for a given normalized intensity and add it to the current final colour (which is 0,0,0)
		Lv += GetEmission(D).ToXYZ();

		// Switch Depending on the shading type
		switch (pScene->m_ShadingType)
		{
			// BRDF Only (Bidirectional Reflectance Distribution Function)
			case 0:
			{
				Lv += UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
				break;
			}
		
			// Phase Function Only
			case 1:
			{
				Lv += 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);
				break;
			}

			// Hybrid (BDRF & Phase Function)
			case 2:
			{
				const float GradMag = GradientMagnitude(Pe) * gIntensityInvRange;

				const float PdfBrdf = (1.0f - __expf(-pScene->m_GradientFactor * GradMag));

				if (RNG.Get1() < PdfBrdf)
  					Lv += UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
				else
 					Lv += 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);

				break;
			}
		}
	}
	else
	{
		if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, INF_MAX), Li, Pl, pLight))
			Lv = Li;
	}

	pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
	//float r = curand_uniform(&localState);
	//pView->m_FrameEstimateXyza.Set(CColorXyza(r, r, r), X, Y);
}

void SingleScattering(CScene* pScene, CScene* pDevScene, CCudaView* pView)
{
	const dim3 KernelBlock(KRNL_SS_BLOCK_W, KRNL_SS_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResX() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResY() / (float)KernelBlock.y));
	
	KrnlSingleScattering<<<KernelGrid, KernelBlock>>>(pDevScene, pView);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Single Scattering");
}	