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
#include "CudaUtilities.h"

#define KRNL_MS_BLOCK_W		32
#define KRNL_MS_BLOCK_H		8
#define KRNL_MS_BLOCK_SIZE	KRNL_MS_BLOCK_W * KRNL_MS_BLOCK_H

DEV bool Terminate(CColorXyz& throughput, CRNG &RNG, float p) {
	if (throughput[0] < 0.01 && throughput[1] < 0.01 && throughput[2] < 0.01) {
		if (RNG.Get1() > p) {
			return true;
		}
		throughput /= 1 - p;
	}
	return false;
}

//KERNEL void KrnlMultipleScattering(CScene* pScene, int* pSeeds)
KERNEL void KrnlMultipleScattering(CScene* pScene, CCudaView* pView)
{
	const int X		= (blockIdx.x * blockDim.x) + threadIdx.x;
	const int Y		= (blockIdx.y * blockDim.y) + threadIdx.y;
	const int PID	= (Y * gFilmWidth) + X;

	if (X >= gFilmWidth || Y >= gFilmHeight || PID >= gFilmNoPixels)
		return;
	
	CRNG RNG(pView->m_RandomSeeds1.GetPtr(X, Y), pView->m_RandomSeeds2.GetPtr(X, Y));

	CColorXyz Lv = SPEC_BLACK, Li = SPEC_BLACK, Tr = SPEC_WHITE;

	CRay Re;
	
	const Vec2f UV = Vec2f(X, Y) + RNG.Get2();

 	pScene->m_Camera.GenerateRay(UV, RNG.Get2(), Re.m_O, Re.m_D);

	Re.m_MinT = 0.0f; 
	Re.m_MaxT = FLT_MAX;

	Vec3f Pe, Pl;
	
	CLight* pLight = NULL;

	for (int i = 0; i < pScene->m_MaxBounces; i++)
	{
		if (SampleDistanceRM(Re, RNG, Pe))
		{
			const float D = GetNormalizedIntensity(Pe);

			Lv += Tr * GetEmission(D).ToXYZ();

			// Switch Depending on the shading type
			switch (pScene->m_ShadingType)
			{
				// BRDF Only (Bidirectional Reflectance Distribution Function)
				case 0:
				{
					Lv += Tr * UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
					break;
				}

				// Phase Function Only
				case 1:
				{
					Lv += Tr * 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);
					break;
				}

				// Hybrid (BDRF & Phase Function)
				case 2:
				{
					const float GradMag = GradientMagnitude(Pe) * gIntensityInvRange;
					const float PdfBrdf = (1.0f - __expf(-pScene->m_GradientFactor * GradMag));
					if (RNG.Get1() < PdfBrdf)
						Lv += Tr * UniformSampleOneLight(pScene, CVolumeShader::Brdf, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, true);
					else
						Lv += Tr * 0.5f * UniformSampleOneLight(pScene, CVolumeShader::Phase, D, Normalize(-Re.m_D), Pe, NormalizedGradient(Pe), RNG, false);

					break;
				}
			}

			// Lets see if we can use the same trick for the direction as for calculating the light
			const float GradMag = GradientMagnitude(Pe) * gIntensityInvRange;
			const float PdfBrdf = (1.0f - __expf(-pScene->m_GradientFactor * GradMag));
			if (RNG.Get1() < PdfBrdf) {
				CVolumeShader Shader(CVolumeShader::Brdf, NormalizedGradient(Pe), Normalize(-Re.m_D), GetDiffuse(D).ToXYZ(), GetSpecular(D).ToXYZ(), 2.5f, GetRoughness(D));
				Vec3f Wi;
				float ShaderPdf = 1.0f;
				CColorXyz F = SPEC_BLACK;
				CLightingSample LS;
				LS.LargeStep(RNG);
				F = Shader.SampleF(Normalize(-Re.m_D), Wi, ShaderPdf, LS.m_BsdfSample);

				Re.m_O = Pe;
				Re.m_D = Wi;
				Re.m_MinT = 0.0f;
				Re.m_MaxT = INF_MAX;

				if (!F.IsBlack() && ShaderPdf > 0)
					Tr *= F * AbsDot(Wi, NormalizedGradient(Pe)) / ShaderPdf;
				else
					break;
			}
			else {
				CVolumeShader Shader(CVolumeShader::Phase, NormalizedGradient(Pe), Normalize(-Re.m_D), GetDiffuse(D).ToXYZ(), GetSpecular(D).ToXYZ(), 2.5f, GetRoughness(D));
				Vec3f Wi;
				float ShaderPdf = 1.0f;
				CColorXyz F = SPEC_BLACK;
				CLightingSample LS;
				LS.LargeStep(RNG);
				F = Shader.SampleF(Normalize(-Re.m_D), Wi, ShaderPdf, LS.m_BsdfSample);

				Re.m_O = Pe;
				Re.m_D = Wi;
				Re.m_MinT = 0.0f;
				Re.m_MaxT = INF_MAX;

				if (!F.IsBlack() && ShaderPdf > 0)
					Tr *= F / ShaderPdf;
			}

			if (Terminate(Tr, RNG, 0.5)) {
				break;
			}
		}
		else
		{
			// If we immediatly miss everything in the volume, use nearestlight to try render lights/background light
			if (i == 0 && NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, INF_MAX), Li, Pl, pLight)) {
					Lv += Tr * Li;
			}
			break;
		}
	}

	__syncthreads();
	
	pView->m_FrameEstimateXyza.Set(CColorXyza(Lv.c[0], Lv.c[1], Lv.c[2]), X, Y);
}

//void MultipleScattering(CScene* pScene, CScene* pDevScene, int* pSeeds)
void MultipleScattering(CScene* pScene, CScene* pDevScene, CCudaView* pView)
{
	//const dim3 KernelBlock(KRNL_MS_BLOCK_W, KRNL_MS_BLOCK_H);
	const dim3 KernelBlock(KRNL_SS_BLOCK_W, KRNL_SS_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResX() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResY() / (float)KernelBlock.y));
	
	KrnlMultipleScattering<<<KernelGrid, KernelBlock>>>(pDevScene, pView);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Multiple Scattering");
}