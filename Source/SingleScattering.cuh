#pragma once

#include "Transport.cuh"

#include <algorithm>
#include <vector>

DEV CColorXyz IncidentLight(CScene* pScene, const Vec3f& Wo, const Vec3f& Pe, const float& D, CCudaRNG& RNG, CSpectral& Spectral)
{
	switch (pScene->m_ShadingType)
	{
		// Do BRDF shading
		case 0:
		{
			return UniformSampleOneLight(pScene, Normalize(Wo), Pe, NormalizedGradient(pScene, Pe), RNG, true, Spectral);
		}
		
		// Do phase function shading
		case 1:
		{
			return 0.5f * UniformSampleOneLight(pScene, Normalize(Wo), Pe, NormalizedGradient(pScene, Pe), RNG, false, Spectral);
			break;
		}

		// Do hybrid shading (BRDF + phase function)
		case 2:
		{
			// Get normalized gradient magnitude
			const float GradMag = GradientMagnitude(pScene, Pe) / pScene->m_GradientMagnitudeRange.GetLength();

			// Get opacity
			const float Tr = GetOpacity(pScene, D)[0];

			// Determine BRDF vs Phase Function scattering
			const float PdfBrdf = 1.0f - __expf(-pScene->m_GradientFactor * GradMag);

			if (RNG.Get1() < PdfBrdf)
			{
  				return UniformSampleOneLight(pScene, Normalize(Wo), Pe, NormalizedGradient(pScene, Pe), RNG, true, Spectral);
			}
			else
			{
				return 0.5f * UniformSampleOneLight(pScene, Normalize(Wo), Pe, NormalizedGradient(pScene, Pe), RNG, false, Spectral);
			}
		}
	}

	return SPEC_BLACK;
}

KERNEL void KrnlSingleScattering(CScene* pScene, unsigned int* pSeeds, CColorXyz* pDevEstFrameXyz)
{
	const int X = (blockIdx.x * blockDim.x) + threadIdx.x;		// Get global y
	const int Y	= (blockIdx.y * blockDim.y) + threadIdx.y;		// Get global x
	
	// Compute sample ID
	const int SID = (Y * (gridDim.x * blockDim.x)) + X;

	// Exit if beyond kernel boundaries
	if (X >= pScene->m_Camera.m_Film.GetWidth() || Y >= pScene->m_Camera.m_Film.GetHeight())
		return;
	
	// Init random number generator
	CCudaRNG RNG(&pSeeds[SID * 2], &pSeeds[SID * 2 + 1]);

	CColorXyz Lv = SPEC_BLACK, Li = SPEC_BLACK;

	CRay Re;

 	// Generate the camera ray
 	pScene->m_Camera.GenerateRay(Vec2f(X, Y), RNG.Get2(), Re.m_O, Re.m_D);

	Re.m_MinT = 0.0f; 
	Re.m_MaxT = FLT_MAX;

	Vec3f Pe, Pl;
	
	CLight* pLight = NULL;

	CSpectral Spectral(false, floorf(RNG.Get1() * 3.0f));

	if (SampleDistanceRM(Re, RNG, Pe, pScene, Spectral))
	{
		if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, (Pe - Re.m_O).Length()), Li, Pl, pLight))
		{
			pDevEstFrameXyz[Y * (int)pScene->m_Camera.m_Film.m_Resolution.GetResX() + X] = Li;
			return;
		}

		// Fetch density
		const float D = Density(pScene, Pe);

		// Add emission
		Lv += GetEmission(pScene, D).ToXYZ();

		// Add incident light
		if (Spectral.m_Enable)
		{
			for (int i = 0; i < 3; i++)
			{
				Spectral.m_Component = i;
				Lv.c[i] += IncidentLight(pScene, -Re.m_D, Pe, D, RNG, Spectral).c[i];
			}
		}
		else
		{
			Lv += IncidentLight(pScene, -Re.m_D, Pe, D, RNG, Spectral);
		}
	}
	else
	{
		if (NearestLight(pScene, CRay(Re.m_O, Re.m_D, 0.0f, INF_MAX), Li, Pl, pLight))
			Lv = Li;
	}

	pDevEstFrameXyz[Y * (int)pScene->m_Camera.m_Film.m_Resolution.GetResX() + X] = Lv;
}

void SingleScattering(CScene* pScene, CScene* pDevScene, unsigned int* pSeeds, CColorXyz* pDevEstFrameXyz)
{
	const dim3 KernelBlock(16, 8);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResX() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.m_Resolution.GetResY() / (float)KernelBlock.y));
	
	KrnlSingleScattering<<<KernelGrid, KernelBlock>>>(pDevScene, pSeeds, pDevEstFrameXyz);
}