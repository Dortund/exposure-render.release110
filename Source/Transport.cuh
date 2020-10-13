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

#include "Shader.cuh"
#include "RayMarching.cuh"

DEV CColorXyz EstimateDirectLight(CScene* pScene, const CVolumeShader::EType& Type, const float& Density, CLight& Light, CLightingSample& LS, const Vec3f& Wo, const Vec3f& Pe, const Vec3f& N, CRNG& RNG)
{
	CColorXyz Ld = SPEC_BLACK, Li = SPEC_BLACK, F = SPEC_BLACK;
	
	CVolumeShader Shader(Type, Pe, N, Wo, GetDiffuse(Density).ToXYZ(), GetSpecular(Density).ToXYZ(), 2.5f, GetRoughness(Density));
	
	CRay Rl; 

	float LightPdf = 1.0f, ShaderPdf = 1.0f;
	
	Vec3f Wi, P, Pl;

	// Take a sample of the light.
	// Also generate a ray 'Rl from a point on the light towards the scattering point
	// Also register the PDF for the light
	// 'LS' ??
 	Li = Light.SampleL(Pe, Rl, LightPdf, LS);
	
	CLight* pLight = NULL;

	// Create vector pointing from scattering point to point on light (direction already normalized)
	Wi = -Rl.m_D; 

	F = Shader.F(Wo, Wi); 

	ShaderPdf = Shader.Pdf(Wo, Wi);
	
	// if the light has a color, PDF's are greater then zero and there is a free enough path between the light point and scattering point
	if (!Li.IsBlack() && ShaderPdf > 0.0f && LightPdf > 0.0f && !FreePathRM(Rl, RNG))
	{
		const float WeightMIS = PowerHeuristic(1.0f, LightPdf, 1.0f, ShaderPdf);
		
		if (Type == CVolumeShader::Brdf)
			Ld += F * Li * AbsDot(Wi, N) * WeightMIS / LightPdf;
		else if (Type == CVolumeShader::Phase)
			Ld += F * Li * WeightMIS / LightPdf;
		else
			Ld += F * Li * WeightMIS / LightPdf;
	}
	
	if (Type == CVolumeShader::OneDirectional || Type == CVolumeShader::LightPathsOctoGradientRejectionSampling)
		F = Shader.SampleFRejection(Wo, Wi, ShaderPdf, LS.m_BsdfSample, RNG);
	else
		F = Shader.SampleF(Wo, Wi, ShaderPdf, LS.m_BsdfSample);

	if (!F.IsBlack() && ShaderPdf > 0.0f)
	{
		if (NearestLight(pScene, CRay(Pe, Wi, 0.0f), Li, Pl, pLight, &LightPdf))
		{
			LightPdf = pLight->Pdf(Pe, Wi);

			if (LightPdf > 0.0f && !Li.IsBlack() && !FreePathRM(CRay(Pl, Normalize(Pe - Pl), 0.0f, (Pe - Pl).Length()), RNG)) 
			{
				const float WeightMIS = PowerHeuristic(1.0f, ShaderPdf, 1.0f, LightPdf);

				if (Type == CVolumeShader::Brdf)
					Ld += F * Li * AbsDot(Wi, N) * WeightMIS / ShaderPdf;
				else if (Type == CVolumeShader::Phase)
					Ld += F * Li * WeightMIS / ShaderPdf;
				else
					Ld += F * Li * WeightMIS / ShaderPdf;
			}
		}
	}

	return Ld;
}

DEV CColorXyz UniformSampleOneLight(CScene* pScene, const CVolumeShader::EType& Type, const float& Density, const Vec3f& Wo, const Vec3f& Pe, const Vec3f& N, CRNG& RNG, const bool& Brdf)
{
	
	const int NumLights = pScene->m_Lighting.m_NoLights;

 	if (NumLights == 0)
 		return SPEC_BLACK;

	CLightingSample LS;
	
	LS.LargeStep(RNG);

	const int WhichLight = (int)floorf(LS.m_LightNum * (float)NumLights);

	CLight& Light = pScene->m_Lighting.m_Lights[WhichLight];
	
	return (float)NumLights * EstimateDirectLight(pScene, Type, Density, Light, LS, Wo, Pe, N, RNG);
}



// Density Property based transport
DEV CColorXyz EstimateDirectLightPropertyBased(CScene* pScene, const CVolumeShader::EType& Type, CLight& Light, CLightingSample& LS, materialProperties &properties, const Vec3f& Wo, const Vec3f& Pe, const Vec3f& N, CRNG& RNG)
{
	CColorXyz Ld = SPEC_BLACK, Li = SPEC_BLACK, F = SPEC_BLACK;

	CVolumeShader Shader(Type, Pe, N, Wo, properties.diffuse.ToXYZ(), properties.specular.ToXYZ(), 2.5f, properties.roughness);

	CRay Rl;

	float LightPdf = 1.0f, ShaderPdf = 1.0f;

	Vec3f Wi, P, Pl;

	// Take a sample of the light.
	// Also generate a ray 'Rl from a point on the light towards the scattering point
	// Also register the PDF for the light
	// 'LS' ??
	Li = Light.SampleL(Pe, Rl, LightPdf, LS);
	
	CLight* pLight = NULL;

	// Create vector pointing from scattering point to point on light (direction already normalized)
	Wi = -Rl.m_D;

	F = Shader.F(Wo, Wi);

	ShaderPdf = Shader.Pdf(Wo, Wi);

	// if the light has a color, PDF's are greater then zero and there is a free enough path between the light point and scattering point
	if (!Li.IsBlack() && ShaderPdf > 0.0f && LightPdf > 0.0f && !FreePathRMPropertyBased(Rl, RNG))
	{
		const float WeightMIS = PowerHeuristic(1.0f, LightPdf, 1.0f, ShaderPdf);

		if (Type == CVolumeShader::Brdf)
			Ld += F * Li * AbsDot(Wi, N) * WeightMIS / LightPdf;
		else if (Type == CVolumeShader::Phase)
			Ld += F * Li * WeightMIS / LightPdf;
		else if (Type == CVolumeShader::LightPaths || Type == CVolumeShader::LightPathsOcto || Type == CVolumeShader::LightPathsOctoGradient)
			Ld += F * Li * WeightMIS / LightPdf;
		else
			Ld += F * Li * WeightMIS / LightPdf;
	}
	
	if (Type == CVolumeShader::OneDirectional || Type == CVolumeShader::LightPathsOctoGradientRejectionSampling)
		F = Shader.SampleFRejection(Wo, Wi, ShaderPdf, LS.m_BsdfSample, RNG);
	else
		F = Shader.SampleF(Wo, Wi, ShaderPdf, LS.m_BsdfSample);

	if (!F.IsBlack() && ShaderPdf > 0.0f)
	{
		if (NearestLight(pScene, CRay(Pe, Wi, 0.0f), Li, Pl, pLight, &LightPdf))
		{
			LightPdf = pLight->Pdf(Pe, Wi);

			if (LightPdf > 0.0f && !Li.IsBlack() && !FreePathRMPropertyBased(CRay(Pl, Normalize(Pe - Pl), 0.0f, (Pe - Pl).Length()), RNG))
			{
				const float WeightMIS = PowerHeuristic(1.0f, ShaderPdf, 1.0f, LightPdf);

				if (Type == CVolumeShader::Brdf)
					Ld += F * Li * AbsDot(Wi, N) * WeightMIS / ShaderPdf;
				else if (Type == CVolumeShader::Phase)
					Ld += F * Li * WeightMIS / ShaderPdf;
				else if (Type == CVolumeShader::LightPaths || Type == CVolumeShader::LightPathsOcto || Type == CVolumeShader::LightPathsOctoGradient)
					Ld += F * Li * WeightMIS / ShaderPdf;
				else
					Ld += F * Li * WeightMIS / ShaderPdf;
			}
		}
	}

	return Ld;
}

DEV CColorXyz UniformSampleOneLightPropertyBased(CScene* pScene, const CVolumeShader::EType& Type, materialProperties &properties, const Vec3f& Wo, const Vec3f& Pe, const Vec3f& N, CRNG& RNG, const bool& Brdf)
{

	const int NumLights = pScene->m_Lighting.m_NoLights;

	if (NumLights == 0)
		return SPEC_BLACK;

	CLightingSample LS;

	LS.LargeStep(RNG);

	const int WhichLight = (int)floorf(LS.m_LightNum * (float)NumLights);

	CLight& Light = pScene->m_Lighting.m_Lights[WhichLight];
	
	return (float)NumLights * EstimateDirectLightPropertyBased(pScene, Type, Light, LS, properties, Wo, Pe, N, RNG);
}