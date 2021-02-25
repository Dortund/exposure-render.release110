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

DEV float GetNormalizedIntensity(const float3& P)
{
	float3 coord = P * gInvAaBbMax * gWorldToTextureTransform - gVoxelTextureOffset;

	if (coord.x < 0 || coord.x >= 1.f ||
		coord.y < 0 || coord.y >= 1.f ||
		coord.z < 0 || coord.z >= 1.f)
		return NAN;

	const float Intensity = ((float)SHRT_MAX * tex3D(gTexDensity, coord.x, coord.y, coord.z));

	return (Intensity - gIntensityMin) * gIntensityInvRange;
}


DEV float GetNormalizedIntensity(const Vec3f& P)
{
	//const float3 c = { P.x, P.y, P.z };
	//return GetNormalizedIntensity(c);
	/*float3 coord = {
		P.x * gInvAaBbMax.x * gWorldToTextureTransform.x,
		P.y * gInvAaBbMax.y * gWorldToTextureTransform.y,
		P.z * gInvAaBbMax.z * gWorldToTextureTransform.z
	};*/

	float3 coord = {
		P.x * gInvAaBbMax.x,// * gWorldToTextureTransform.x - gVoxelTextureOffset.x,
		P.y * gInvAaBbMax.y,// * gWorldToTextureTransform.y - gVoxelTextureOffset.y,
		P.z * gInvAaBbMax.z,// * gWorldToTextureTransform.z - gVoxelTextureOffset.z
	};

	/*if (coord.x < 0 || coord.x > gWorldToTextureTransform.x ||
		coord.y < 0 || coord.y > gWorldToTextureTransform.y ||
		coord.z < 0 || coord.z > gWorldToTextureTransform.z)
		return NAN;*/
	/*if (coord.x < -gVoxelTextureOffset.x / 2 || coord.x > 1.f - gVoxelTextureOffset.x / 2 ||
		coord.y < -gVoxelTextureOffset.y / 2 || coord.y > 1.f - gVoxelTextureOffset.y / 2 ||
		coord.z < -gVoxelTextureOffset.z / 2 || coord.z > 1.f - gVoxelTextureOffset.z / 2)
		return NAN;*/

	const float Intensity = ((float)SHRT_MAX * tex3D(gTexDensity, coord.x, coord.y, coord.z));

	return (Intensity - gIntensityMin) * gIntensityInvRange;
}

DEV float GetOpacity(const float& NormalizedIntensity)
{
	if (isnan(NormalizedIntensity))
		return 0;

	return tex1D(gTexOpacity, NormalizedIntensity);
}

DEV CColorRgbHdr GetDiffuse(const float& NormalizedIntensity)
{
	if (isnan(NormalizedIntensity))
		return CColorRgbHdr(0,0,0);

	float4 Diffuse = tex1D(gTexDiffuse, NormalizedIntensity);
	return CColorRgbHdr(Diffuse.x, Diffuse.y, Diffuse.z);
}

DEV CColorRgbHdr GetSpecular(const float& NormalizedIntensity)
{
	if (isnan(NormalizedIntensity))
		return CColorRgbHdr(0);

	float4 Specular = tex1D(gTexSpecular, NormalizedIntensity);
	return CColorRgbHdr(Specular.x, Specular.y, Specular.z);
}

DEV float GetRoughness(const float& NormalizedIntensity)
{
	if (isnan(NormalizedIntensity))
		return 0;

	return tex1D(gTexRoughness, NormalizedIntensity);
}

DEV CColorRgbHdr GetEmission(const float& NormalizedIntensity)
{
	if (isnan(NormalizedIntensity))
		return CColorRgbHdr(0);

	float4 Emission = tex1D(gTexEmission, NormalizedIntensity);
	return CColorRgbHdr(Emission.x, Emission.y, Emission.z);
}

DEV inline Vec3f NormalizedGradient(const Vec3f& P)
{
	Vec3f Gradient;

	Gradient.x = (GetOpacity(GetNormalizedIntensity(P + ToVec3f(gVoxelSizeWorldX))) - GetOpacity(GetNormalizedIntensity(P - ToVec3f(gVoxelSizeWorldX)))) * gInvGradientDelta.x;
	Gradient.y = (GetOpacity(GetNormalizedIntensity(P + ToVec3f(gVoxelSizeWorldY))) - GetOpacity(GetNormalizedIntensity(P - ToVec3f(gVoxelSizeWorldY)))) * gInvGradientDelta.y;
	Gradient.z = (GetOpacity(GetNormalizedIntensity(P + ToVec3f(gVoxelSizeWorldZ))) - GetOpacity(GetNormalizedIntensity(P - ToVec3f(gVoxelSizeWorldZ)))) * gInvGradientDelta.z;

	return Normalize(Gradient);
}

DEV float GradientMagnitude(const Vec3f& P)
{
	return ((float)SHRT_MAX * tex3D(gTexGradientMagnitude, P.x * gInvAaBbMax.x, P.y * gInvAaBbMax.y, P.z * gInvAaBbMax.z));
}

DEV int GetLightPathValue(const float3 P) {
	float3 coord = P * gInvAaBbMax * gWorldToTextureTransform - gVoxelTextureOffset;

	/*if (coord.x < 0 || coord.x >= 1.f ||
		coord.y < 0 || coord.y >= 1.f ||
		coord.z < 0 || coord.z >= 1.f)
		return NAN;*/

	return tex3D(gTexLightPaths, coord.x, coord.y, coord.z);
}

/*
DEV inline float biLinear(float tx, float ty, float c00, float c10, float c01, float c11) {
	return (1 - tx) * (1 - ty) * c00 +
		tx * (1 - ty) * c10 +
		(1 - tx) * ty * c01 +
		tx * ty * c11;
}

DEV inline CColorRgbHdr biLinear(float tx, float ty,const CColorRgbHdr& c00, const CColorRgbHdr& c10, const CColorRgbHdr& c01, const CColorRgbHdr& c11) {
	CColorRgbHdr  a = Lerp(tx, c00, c10);
	CColorRgbHdr  b = Lerp(tx, c01, c11);
	return Lerp(ty, a, b);
}
*/
struct materialProperties {
	float opacity;
	float roughness;
	CColorRgbHdr diffuse;
	CColorRgbHdr specular;
	CColorRgbHdr emission;
};

DEV inline void blendMaterialProperties(materialProperties &result, materialProperties &a, materialProperties &b, float amount)
{
	//float c = 100;
	//float powAmount = pow(amount, ((c * gDensityScale * b.opacity + 1) / (c * gDensityScale * a.opacity + 1)));
	float powAmount = amount;
	
	//result.opacity = lerp(a.opacity, b.opacity, powAmount);
	

	/*if (gFilmWidth == 500) {
		result.opacity = lerp(a.opacity, b.opacity, amount);
		result.roughness = lerp(a.roughness * a.opacity, b.roughness * b.opacity, powAmount);
		result.diffuse = Lerp(powAmount, a.diffuse * a.opacity, b.diffuse * b.opacity);
		result.specular = Lerp(powAmount, a.specular * a.opacity, b.specular * b.opacity);
		result.emission = Lerp(powAmount, a.emission * a.opacity, b.emission * b.opacity);
	}
	else if (gFilmWidth == 501) {
		result.opacity = lerp(a.opacity, b.opacity, amount);
		result.roughness = lerp(a.roughness, b.roughness, powAmount);
		result.diffuse = Lerp(powAmount, a.diffuse, b.diffuse);
		result.specular = Lerp(powAmount, a.specular, b.specular);
		result.emission = Lerp(powAmount, a.emission, b.emission);
	}
	else if (gFilmWidth == 502) {
		//float opMin = fmin(a.opacity, b.opacity);
		//float opMax = fmax(a.opacity, b.opacity);
		*/
		//powAmount = (a.opacity - opMin) / (opMax - opMin);
		if (a.opacity > 0 || b.opacity > 0) {
			float opAmount = b.opacity / (a.opacity + b.opacity);
			powAmount = lerp(powAmount, opAmount, abs(opAmount - 0.5) * 2);
		}

		result.opacity = lerp(a.opacity, b.opacity, amount);
		result.roughness = lerp(a.roughness, b.roughness, powAmount);
		result.diffuse = Lerp(powAmount, a.diffuse, b.diffuse);
		result.specular = Lerp(powAmount, a.specular, b.specular);
		result.emission = Lerp(powAmount, a.emission, b.emission);
	/*}
	else {
		if (a.opacity == 0)
			powAmount = 1;
		if (b.opacity == 0)
			powAmount = 0;

		result.opacity = lerp(a.opacity, b.opacity, amount);
		result.roughness = lerp(a.roughness, b.roughness, powAmount);
		result.diffuse = Lerp(powAmount, a.diffuse, b.diffuse);
		result.specular = Lerp(powAmount, a.specular, b.specular);
		result.emission = Lerp(powAmount, a.emission, b.emission);
	}*/
}

DEV inline void sampleRawProperties(materialProperties &properties, const float3& P)
{
	float intensity = GetNormalizedIntensity(P);
	properties.opacity = GetOpacity(intensity);
	properties.roughness = GetRoughness(intensity);
	properties.diffuse = GetDiffuse(intensity);
	properties.specular = GetSpecular(intensity);
	properties.emission = GetEmission(intensity);
}

DEV void sampleProperties(materialProperties &properties, float3 *opacityGradient, float3 *fractions, const Vec3f& P)
{
	/*
	
	[I]=G________F  +y
	    /      /|    ^
	   /      / |    |
	 H/_____E/  |
	  |  C  |  /B       +z
	  |     | /        /
	  |_____|/       -z
	 D     A=(0,0,0)
	    +x<->-x
	*/

	float3 coord;
	coord.x = P.x;
	coord.y = P.y;
	coord.z = P.z;

	//A is reference
	float3 mins = floor(coord / gVoxelSizeWorld) * gVoxelSizeWorld;
	float3 fraction = clamp((coord - mins) / gVoxelSizeWorld, 0, 1);

	if (fractions) {
		fractions->x = fraction.x;
		fractions->y = fraction.y;
		fractions->z = fraction.z;
	}

	materialProperties A, B, C, D, E, F, G, H;
	sampleRawProperties(A, mins);
	sampleRawProperties(B, mins + gVoxelSizeWorldZ);
	sampleRawProperties(C, mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
	sampleRawProperties(D, mins + gVoxelSizeWorldX);

	sampleRawProperties(E, mins + gVoxelSizeWorldY);
	sampleRawProperties(F, mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
	sampleRawProperties(G, mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
	sampleRawProperties(H, mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

	materialProperties AD, BC, EH, FG;
	blendMaterialProperties(AD, A, D, fraction.x);
	blendMaterialProperties(BC, B, C, fraction.x);
	blendMaterialProperties(EH, E, H, fraction.x);
	blendMaterialProperties(FG, F, G, fraction.x);

	materialProperties ADBC, EHFG;
	blendMaterialProperties(ADBC, AD, BC, fraction.z);
	blendMaterialProperties(EHFG, EH, FG, fraction.z);

	blendMaterialProperties(properties, ADBC, EHFG, fraction.y);
	
	if (opacityGradient)
	{
		float AB, EF, DC, HG, ABEF, DCHG;
		AB = lerp(A.opacity, B.opacity, fraction.z);
		EF = lerp(E.opacity, F.opacity, fraction.z);
		DC = lerp(D.opacity, C.opacity, fraction.z);
		HG = lerp(H.opacity, G.opacity, fraction.z);
		ABEF = lerp(AB, EF, fraction.y);
		DCHG = lerp(DC, HG, fraction.y);
		opacityGradient->x = (DCHG - ABEF) / gSpacings.x;

		opacityGradient->y = (EHFG.opacity - ADBC.opacity) / gSpacings.y;

		float ADEH, BCFG;
		//AD = lerp(A.opacity, D.opacity, fraction.x);
		//EH = lerp(E.opacity, H.opacity, fraction.x);
		//BC = lerp(B.opacity, C.opacity, fraction.x);
		//FG = lerp(F.opacity, G.opacity, fraction.x);
		ADEH = lerp(AD.opacity, EH.opacity, fraction.y);
		BCFG = lerp(BC.opacity, FG.opacity, fraction.y);
		opacityGradient->z = (BCFG - ADEH) / gSpacings.z;
	}
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

	 //return A + ((Ax - A) / max((float)N, 1.0f));
	// N = 0: A + ((AX - A) / 1) == ? + ((x - ?) / 1) == x
	// N = 1: A + ((AX - A) / 1) == x + ((x2 - x) / 1) == x + x2 - x == x2
	// N = 2: A + ((AX - A) / 2) == x2 + ((x3 - x2) / 2) == x2 + (x3/2) - (x2/2) == (x2 + x3) / 2

	return A + ((Ax - A) / (N+1));
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

HOD Vec3i Index1To3(int i, int sizeX, int sizeY, int sizeZ) {
	Vec3i result;
	int XY = (sizeX * sizeY);
	result.z = i / XY;
	int remainder = i - XY * result.z;
	result.y = remainder / sizeX;
	result.x = remainder % sizeX;
	return result;
}