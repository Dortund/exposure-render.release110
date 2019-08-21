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
	float3 coord = {
		P.x * gInvTextureSize.x,
		P.y * gInvTextureSize.y,
		P.z * gInvTextureSize.z
	};

	if (coord.x < 0 || coord.x > gInvTextureSize.x ||
		coord.y < 0 || coord.y > gInvTextureSize.y ||
		coord.z < 0 || coord.z > gInvTextureSize.z)
		return NAN;

	const float Intensity = ((float)SHRT_MAX * tex3D(gTexDensity, coord.x, coord.y, coord.z));

	return (Intensity - gIntensityMin) * gIntensityInvRange;
}
DEV float GetNormalizedIntensity(const float3& P)
{
	float3 coord = P * gInvTextureSize;

	if (coord.x < 0 || coord.x > gInvTextureSize.x ||
		coord.y < 0 || coord.y > gInvTextureSize.y ||
		coord.z < 0 || coord.z > gInvTextureSize.z)
		return NAN;

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

struct materialProperties {
	float opacity;
	float roughness;
	CColorRgbHdr diffuse;
	CColorRgbHdr specular;
	CColorRgbHdr emission;
};

DEV void blendMaterialProperties(materialProperties &result, materialProperties &a, materialProperties &b, float amount)
{
	result.opacity = lerp(a.opacity, b.opacity, amount);
	result.roughness = lerp(a.roughness, b.roughness, amount);
	result.diffuse = Lerp(amount, a.diffuse, b.diffuse);
	result.specular = Lerp(amount, a.specular, b.specular);
	result.emission = Lerp(amount, a.emission, b.emission);
}

DEV void sampleRawProperties(materialProperties &properties, const float3& P)
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
	   F________G
	   /      /|
	  /      / |
    E_______/H |
	 |   B |  /C
	 |     | /
	 |_____|/
	 A     D
	 bottom
	*/

	float3 coord;
	coord.x = P.x;
	coord.y = P.y;
	coord.z = P.z;

	//coord = (coord - 0.5 * gGradientDelta) / (1 - gGradientDelta);

	//A is reference
	float3 mins = floor(coord / gGradientDelta) * gGradientDelta;
	float3 fraction = (coord - mins) / gGradientDelta;

	if (fractions) {
		fractions->x = fraction.x;
		fractions->y = fraction.y;
		fractions->z = fraction.z;
	}

	materialProperties A, B, C, D, E, F, G, H;
	sampleRawProperties(A, mins);
	sampleRawProperties(B, mins + gGradientDeltaZ);
	sampleRawProperties(C, mins + gGradientDeltaZ + gGradientDeltaX);
	sampleRawProperties(D, mins + gGradientDeltaX);

	sampleRawProperties(E, mins + gGradientDeltaY);
	sampleRawProperties(F, mins + gGradientDeltaZ + gGradientDeltaY);
	sampleRawProperties(G, mins + gGradientDeltaZ + gGradientDeltaX + gGradientDeltaY);
	sampleRawProperties(H, mins + gGradientDeltaX + gGradientDeltaY);

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
		opacityGradient->y = EHFG.opacity - ADBC.opacity;
		
		float AB, EF, DC, HG, ABEF, DCHG;
		AB = lerp(A.opacity, B.opacity, fraction.z);
		EF = lerp(E.opacity, F.opacity, fraction.z);
		DC = lerp(D.opacity, C.opacity, fraction.z);
		HG = lerp(H.opacity, G.opacity, fraction.z);
		ABEF = lerp(AB, EF, fraction.y);
		DCHG = lerp(DC, HG, fraction.y);
		opacityGradient->x = DCHG - ABEF;

		float ADEH, BCFG;
		//AD = lerp(A.opacity, D.opacity, fraction.x);
		//EH = lerp(E.opacity, H.opacity, fraction.x);
		//BC = lerp(B.opacity, C.opacity, fraction.x);
		//FG = lerp(F.opacity, G.opacity, fraction.x);
		ADEH = lerp(AD.opacity, EH.opacity, fraction.y);
		BCFG = lerp(BC.opacity, FG.opacity, fraction.y);
		opacityGradient->z = BCFG - ADEH;
	}
}



DEV float GetOpacityProperty(const Vec3f& P) {
	float gxi, gyi, gzi, tx, ty, tz;
	gxi = int(P.x / gGradientDeltaX.x) * gGradientDeltaX.x; tx = (P.x - gxi) / gGradientDeltaX.x;
	gyi = int(P.y / gGradientDeltaY.y) * gGradientDeltaY.y; ty = (P.y - gyi) / gGradientDeltaY.y;
	gzi = int(P.z / gGradientDeltaZ.z) * gGradientDeltaZ.z; tz = (P.z - gzi) / gGradientDeltaZ.z;
	const float c000 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi)));
	const float c100 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi)));
	const float c010 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi)));
	const float c110 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi)));
	const float c001 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi + gGradientDeltaZ.z)));
	const float c101 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi + gGradientDeltaZ.z)));
	const float c011 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));
	const float c111 = GetOpacity(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));

	float e = biLinear(tx, ty, c000, c100, c010, c110);
	float f = biLinear(tx, ty, c001, c101, c011, c111);
	//return lerp(tz, f, e);
	return e * (1 - tz) + f * tz;
}

DEV inline Vec3f GetNormalizedOpacityGradientProperty(const Vec3f& P)
{
	Vec3f Gradient;

	Gradient.x = (GetOpacityProperty(P + ToVec3f(gGradientDeltaX)) - GetOpacityProperty(P - ToVec3f(gGradientDeltaX))) * gInvGradientDelta;
	Gradient.y = (GetOpacityProperty(P + ToVec3f(gGradientDeltaY)) - GetOpacityProperty(P - ToVec3f(gGradientDeltaY))) * gInvGradientDelta;
	Gradient.z = (GetOpacityProperty(P + ToVec3f(gGradientDeltaZ)) - GetOpacityProperty(P - ToVec3f(gGradientDeltaZ))) * gInvGradientDelta;

	return Normalize(Gradient);
}

DEV inline float GetOpacityGradientMagnitudeProperty(const Vec3f& P)
{
	Vec3f Gradient;

	Gradient.x = (GetOpacityProperty(P + ToVec3f(gGradientDeltaX)) - GetOpacityProperty(P - ToVec3f(gGradientDeltaX))) * gInvGradientDelta;
	Gradient.y = (GetOpacityProperty(P + ToVec3f(gGradientDeltaY)) - GetOpacityProperty(P - ToVec3f(gGradientDeltaY))) * gInvGradientDelta;
	Gradient.z = (GetOpacityProperty(P + ToVec3f(gGradientDeltaZ)) - GetOpacityProperty(P - ToVec3f(gGradientDeltaZ))) * gInvGradientDelta;

	return Gradient.Length();
}

DEV CColorRgbHdr GetDiffuseProperty(const Vec3f& P) {
	float gxi, gyi, gzi, tx, ty, tz;
	gxi = int(P.x / gGradientDeltaX.x) * gGradientDeltaX.x; tx = (P.x - gxi) / gGradientDeltaX.x;
	gyi = int(P.y / gGradientDeltaY.y) * gGradientDeltaY.y; ty = (P.y - gyi) / gGradientDeltaY.y;
	gzi = int(P.z / gGradientDeltaZ.z) * gGradientDeltaZ.z; tz = (P.z - gzi) / gGradientDeltaZ.z;
	const CColorRgbHdr c000 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi)));
	const CColorRgbHdr c100 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi)));
	const CColorRgbHdr c010 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi)));
	const CColorRgbHdr c110 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi)));
	const CColorRgbHdr c001 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c101 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c011 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c111 = GetDiffuse(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));

	CColorRgbHdr e = biLinear(tx, ty, c000, c100, c010, c110);
	CColorRgbHdr f = biLinear(tx, ty, c001, c101, c011, c111);
	return e * (1 - tz) + f * tz;
}

DEV CColorRgbHdr GetSpecularProperty(const Vec3f& P) {
	float gxi, gyi, gzi, tx, ty, tz;
	gxi = int(P.x / gGradientDeltaX.x) * gGradientDeltaX.x; tx = (P.x - gxi) / gGradientDeltaX.x;
	gyi = int(P.y / gGradientDeltaY.y) * gGradientDeltaY.y; ty = (P.y - gyi) / gGradientDeltaY.y;
	gzi = int(P.z / gGradientDeltaZ.z) * gGradientDeltaZ.z; tz = (P.z - gzi) / gGradientDeltaZ.z;
	const CColorRgbHdr c000 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi)));
	const CColorRgbHdr c100 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi)));
	const CColorRgbHdr c010 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi)));
	const CColorRgbHdr c110 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi)));
	const CColorRgbHdr c001 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c101 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c011 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c111 = GetSpecular(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));

	CColorRgbHdr e = biLinear(tx, ty, c000, c100, c010, c110);
	CColorRgbHdr f = biLinear(tx, ty, c001, c101, c011, c111);
	return e * (1 - tz) + f * tz;
}

DEV float GetRoughnessProperty(const Vec3f& P) {
	float gxi, gyi, gzi, tx, ty, tz;
	gxi = int(P.x / gGradientDeltaX.x) * gGradientDeltaX.x; tx = (P.x - gxi) / gGradientDeltaX.x;
	gyi = int(P.y / gGradientDeltaY.y) * gGradientDeltaY.y; ty = (P.y - gyi) / gGradientDeltaY.y;
	gzi = int(P.z / gGradientDeltaZ.z) * gGradientDeltaZ.z; tz = (P.z - gzi) / gGradientDeltaZ.z;
	const float c000 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi)));
	const float c100 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi)));
	const float c010 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi)));
	const float c110 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi)));
	const float c001 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi + gGradientDeltaZ.z)));
	const float c101 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi + gGradientDeltaZ.z)));
	const float c011 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));
	const float c111 = GetRoughness(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));

	float e = biLinear(tx, ty, c000, c100, c010, c110);
	float f = biLinear(tx, ty, c001, c101, c011, c111);
	return e * (1 - tz) + f * tz;
}

DEV CColorRgbHdr GetEmmisionProperty(const Vec3f& P) {
	float gxi, gyi, gzi, tx, ty, tz;
	gxi = int(P.x / gGradientDeltaX.x) * gGradientDeltaX.x; tx = (P.x - gxi) / gGradientDeltaX.x;
	gyi = int(P.y / gGradientDeltaY.y) * gGradientDeltaY.y; ty = (P.y - gyi) / gGradientDeltaY.y;
	gzi = int(P.z / gGradientDeltaZ.z) * gGradientDeltaZ.z; tz = (P.z - gzi) / gGradientDeltaZ.z;
	const CColorRgbHdr c000 = GetEmission(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi)));
	const CColorRgbHdr c100 = GetEmission(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi)));
	const CColorRgbHdr c010 = GetEmission(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi)));
	const CColorRgbHdr c110 = GetEmission(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi)));
	const CColorRgbHdr c001 = GetEmission(GetNormalizedIntensity(Vec3f(gxi, gyi, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c101 = GetEmission(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c011 = GetEmission(GetNormalizedIntensity(Vec3f(gxi, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));
	const CColorRgbHdr c111 = GetEmission(GetNormalizedIntensity(Vec3f(gxi + gGradientDeltaX.x, gyi + gGradientDeltaY.y, gzi + gGradientDeltaZ.z)));

	CColorRgbHdr e = biLinear(tx, ty, c000, c100, c010, c110);
	CColorRgbHdr f = biLinear(tx, ty, c001, c101, c011, c111);
	return e * (1 - tz) + f * tz;
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