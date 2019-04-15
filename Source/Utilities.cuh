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

#include <cuda_runtime.h>

DEV inline Vec3f ToVec3f(const float3& V)
{
	return Vec3f(V.x, V.y, V.z);
}

DEV float GetNormalizedIntensity(const Vec3f& P)
{
	const float Intensity = ((float)SHRT_MAX * tex3D(gTexDensity, P.x * gInvAaBbMax.x, P.y * gInvAaBbMax.y, P.z * gInvAaBbMax.z));

	return (Intensity - gIntensityMin) * gIntensityInvRange;
}

DEV float GetOpacity(const float& NormalizedIntensity)
{
	return tex1D(gTexOpacity, NormalizedIntensity);
}

DEV CColorRgbHdr GetDiffuse(const float& NormalizedIntensity)
{
	float4 Diffuse = tex1D(gTexDiffuse, NormalizedIntensity);
	return CColorRgbHdr(Diffuse.x, Diffuse.y, Diffuse.z);
}

DEV CColorRgbHdr GetSpecular(const float& NormalizedIntensity)
{
	float4 Specular = tex1D(gTexSpecular, NormalizedIntensity);
	return CColorRgbHdr(Specular.x, Specular.y, Specular.z);
}

DEV float GetRoughness(const float& NormalizedIntensity)
{
	return tex1D(gTexRoughness, NormalizedIntensity);
}

DEV CColorRgbHdr GetEmission(const float& NormalizedIntensity)
{
	float4 Emission = tex1D(gTexEmission, NormalizedIntensity);
	return CColorRgbHdr(Emission.x, Emission.y, Emission.z);
}

DEV inline Vec3f NormalizedGradient(const Vec3f& P)
{
	Vec3f Gradient;

	Gradient.x = (GetOpacity(GetNormalizedIntensity(P + ToVec3f(gGradientDeltaX))) - GetOpacity(GetNormalizedIntensity(P - ToVec3f(gGradientDeltaX)))) * gInvGradientDelta;
	Gradient.y = (GetOpacity(GetNormalizedIntensity(P + ToVec3f(gGradientDeltaY))) - GetOpacity(GetNormalizedIntensity(P - ToVec3f(gGradientDeltaY)))) * gInvGradientDelta;
	Gradient.z = (GetOpacity(GetNormalizedIntensity(P + ToVec3f(gGradientDeltaZ))) - GetOpacity(GetNormalizedIntensity(P - ToVec3f(gGradientDeltaZ)))) * gInvGradientDelta;

	return Normalize(Gradient);
}

DEV float GradientMagnitude(const Vec3f& P)
{
	return ((float)SHRT_MAX * tex3D(gTexGradientMagnitude, P.x * gInvAaBbMax.x, P.y * gInvAaBbMax.y, P.z * gInvAaBbMax.z));
}

// TODO is called NearestLight, but seems to take the last of all lights in the list it had a hit with
DEV bool NearestLight(CScene* pScene, CRay R, CColorXyz& LightColor, Vec3f& Pl, CLight*& pLight, float* pPdf = NULL)
{
	bool Hit = false;
	
	float T = 0.0f;

	CRay RayCopy = R;

	float Pdf = 0.0f;

	// Loop over all the lights in the scene
	for (int i = 0; i < pScene->m_Lighting.m_NoLights; i++)
	{
		// Tries to intersect a light with the ray
		// T = distance, gets filled during the call
		// LightColor = light color, gets filled during the call
		// Pdf = probability density function, a float so just chance I guess..., get filled during the call
		if (pScene->m_Lighting.m_Lights[i].Intersect(RayCopy, T, LightColor, NULL, &Pdf))
		{
			Pl		= R(T);
			pLight	= &pScene->m_Lighting.m_Lights[i];
			Hit		= true;
		}
	}
	
	if (pPdf)
		*pPdf = Pdf;

	return Hit;
}

/// <summary>
/// Checks of the ray intersects with the boundingbox of the volume and fills minT and maxT with the correct values
/// </summary>
DEV bool IntersectBox(const CRay& R, float* pNearT, float* pFarT)
{
	const Vec3f InvR		= Vec3f(1.0f, 1.0f, 1.0f) / R.m_D;
	const Vec3f BottomT		= InvR * (Vec3f(gAaBbMin.x, gAaBbMin.y, gAaBbMin.z) - R.m_O);
	const Vec3f TopT		= InvR * (Vec3f(gAaBbMax.x, gAaBbMax.y, gAaBbMax.z) - R.m_O);
	const Vec3f MinT		= MinVec3f(TopT, BottomT);
	const Vec3f MaxT		= MaxVec3f(TopT, BottomT);
	const float LargestMinT = fmaxf(fmaxf(MinT.x, MinT.y), fmaxf(MinT.x, MinT.z));
	const float LargestMaxT = fminf(fminf(MaxT.x, MaxT.y), fminf(MaxT.x, MaxT.z));
	
	*pNearT = LargestMinT;
	*pFarT	= LargestMaxT;
	
	return LargestMaxT > LargestMinT;
}

DEV CColorXyza CumulativeMovingAverage(const CColorXyza& A, const CColorXyza& Ax, const int& N)
{
//	if (gNoIterations == 0)
//		return CColorXyza(0.0f);

	 return A + ((Ax - A) / max((float)N, 1.0f));
}

HOD int Index3To1(int x, int y, int z, int sizeX, int sizeY, int sizeZ) {
	return z * sizeX * sizeY + y * sizeX + x;
}

HOD int Index3To1(int x, int y, int z, CResolution3D resolution) {
	return Index3To1(x, y, z, resolution.GetResX(), resolution.GetResY(), resolution.GetResZ());
}

HOD int Index3To1(Vec3i p, CResolution3D resolution) {
	return Index3To1(p.x, p.y, p.z, resolution.GetResX(), resolution.GetResY(), resolution.GetResZ());
}

HOD Vec3i Inedex1To3(int i, int sizeX, int sizeY, int sizeZ) {
	Vec3i result;
	int XY = (sizeX * sizeY);
	result.z = i / XY;
	int remainder = i - XY * result.z;
	result.y = remainder / sizeX;
	result.x = remainder % sizeX;
	return result;
}