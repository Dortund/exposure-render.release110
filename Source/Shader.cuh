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

#include "Geometry.h"

#include "MonteCarlo.cuh"
#include "Sample.cuh"

class CLambertian
{
public:
	HOD CLambertian(const CColorXyz& Kd)
	{
		m_Kd = Kd;
	}

	HOD ~CLambertian(void)
	{
	}

	HOD CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_PI_F;
	}

	HOD CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec2f& U)
	{
		Wi = CosineWeightedHemisphere(U);

		if (Wo.z < 0.0f)
			Wi.z *= -1.0f;

		Pdf = this->Pdf(Wo, Wi);

		return this->F(Wo, Wi);
	}

	HOD float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		return SameHemisphere(Wo, Wi) ? AbsCosTheta(Wi) * INV_PI_F : 0.0f;
	}

	CColorXyz	m_Kd;
};

HOD inline CColorXyz FrDiel(float cosi, float cost, const CColorXyz &etai, const CColorXyz &etat)
{
	CColorXyz Rparl = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
	CColorXyz Rperp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
	return (Rparl*Rparl + Rperp*Rperp) / 2.f;
}

class CFresnel
{
public:
	HOD CFresnel(float ei, float et) :
	  eta_i(ei),
		  eta_t(et)
	  {
	  }

	  HOD  ~CFresnel(void)
	  {
	  }

	  HOD CColorXyz Evaluate(float cosi)
	  {
		  // Compute Fresnel reflectance for dielectric
		  cosi = Clamp(cosi, -1.0f, 1.0f);

		  // Compute indices of refraction for dielectric
		  bool entering = cosi > 0.0f;
		  float ei = eta_i, et = eta_t;

		  if (!entering)
			  swap(ei, et);

		  // Compute _sint_ using Snell's law
		  float sint = ei/et * sqrtf(max(0.f, 1.f - cosi*cosi));

		  if (sint >= 1.0f)
		  {
			  // Handle total internal reflection
			  return 1.0f;
		  }
		  else
		  {
			  float cost = sqrtf(max(0.f, 1.0f - sint * sint));
			  return FrDiel(fabsf(cosi), cost, ei, et);
		  }
	  }

	  float eta_i, eta_t;
};

class CBlinn
{
public:
	HOD CBlinn(const float& Exponent) :
	  m_Exponent(Exponent)
	  {
	  }

	  HOD ~CBlinn(void)
	  {
	  }

	  HOD void SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec2f& U)
	  {
		  // Compute sampled half-angle vector $\wh$ for Blinn distribution
		  float costheta = powf(U.x, 1.f / (m_Exponent+1));
		  float sintheta = sqrtf(max(0.f, 1.f - costheta*costheta));
		  float phi = U.y * 2.f * PI_F;

		  Vec3f wh = SphericalDirection(sintheta, costheta, phi);

		  if (!SameHemisphere(Wo, wh))
			  wh = -wh;

		  // Compute incident direction by reflecting about $\wh$
		  Wi = -Wo + 2.f * Dot(Wo, wh) * wh;

		  // Compute PDF for $\wi$ from Blinn distribution
		  float blinn_pdf = ((m_Exponent + 1.f) * powf(costheta, m_Exponent)) / (2.f * PI_F * 4.f * Dot(Wo, wh));

		  if (Dot(Wo, wh) <= 0.f)
			  blinn_pdf = 0.f;

		  Pdf = blinn_pdf;
	  }

	  HOD float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	  {
		  Vec3f wh = Normalize(Wo + Wi);

		  float costheta = AbsCosTheta(wh);
		  // Compute PDF for $\wi$ from Blinn distribution
		  float blinn_pdf = ((m_Exponent + 1.f) * powf(costheta, m_Exponent)) / (2.f * PI_F * 4.f * Dot(Wo, wh));

		  if (Dot(Wo, wh) <= 0.0f)
			  blinn_pdf = 0.0f;

		  return blinn_pdf;
	  }

	  HOD float D(const Vec3f& wh)
	  {
		  float costhetah = AbsCosTheta(wh);
		  return (m_Exponent+2) * INV_TWO_PI_F * powf(costhetah, m_Exponent);
	  }

	  float	m_Exponent;
};

class CMicrofacet
{
public:
	HOD CMicrofacet(const CColorXyz& Reflectance, const float& Ior, const float& Exponent) :
	  m_R(Reflectance),
		  m_Fresnel(Ior, 1.0f),
		  m_Blinn(Exponent)
	  {
	  }

	  HOD ~CMicrofacet(void)
	  {
	  }

	  HOD CColorXyz F(const Vec3f& wo, const Vec3f& wi)
	  {
		  float cosThetaO = AbsCosTheta(wo);
		  float cosThetaI = AbsCosTheta(wi);

		  if (cosThetaI == 0.f || cosThetaO == 0.f)
			  return SPEC_BLACK;

		  Vec3f wh = wi + wo;

		  if (wh.x == 0. && wh.y == 0. && wh.z == 0.)
			  return SPEC_BLACK;

		  wh = Normalize(wh);
		  float cosThetaH = Dot(wi, wh);

		  CColorXyz F = SPEC_WHITE;//m_Fresnel.Evaluate(cosThetaH);

		  //probably Torrance–Sparrow BRDF
		  return m_R * m_Blinn.D(wh) * G(wo, wi, wh) * F / (4.f * cosThetaI * cosThetaO);
	  }

	  HOD CColorXyz SampleF(const Vec3f& wo, Vec3f& wi, float& Pdf, const Vec2f& U)
	  {
		  m_Blinn.SampleF(wo, wi, Pdf, U);

		  if (!SameHemisphere(wo, wi))
			  return SPEC_BLACK;

		  return this->F(wo, wi);
	  }

	  HOD float Pdf(const Vec3f& wo, const Vec3f& wi)
	  {
		  if (!SameHemisphere(wo, wi))
			  return 0.0f;

		  return m_Blinn.Pdf(wo, wi);
	  }

	  HOD float G(const Vec3f& wo, const Vec3f& wi, const Vec3f& wh)
	  {
		  float NdotWh = AbsCosTheta(wh);
		  float NdotWo = AbsCosTheta(wo);
		  float NdotWi = AbsCosTheta(wi);
		  float WOdotWh = AbsDot(wo, wh);

		  return min(1.f, min((2.f * NdotWh * NdotWo / WOdotWh), (2.f * NdotWh * NdotWi / WOdotWh)));
	  }

	  CColorXyz		m_R;
	  CFresnel		m_Fresnel;
	  CBlinn		m_Blinn;

};

class CLightPaths
{
public:
	DEV CLightPaths(const Vec3f& Pe, const CColorXyz& Kd) :
		m_Kd(Kd),
		m_Pe(Pe)
	{
	}

	DEV ~CLightPaths(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec3f& U)
	{
		/*float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = (1 + GetLightPathValue(mins));
		B = (1 + GetLightPathValue(mins + gVoxelSizeWorldZ));
		C = (1 + GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX));
		D = (1 + GetLightPathValue(mins + gVoxelSizeWorldX));

		E = (1 + GetLightPathValue(mins + gVoxelSizeWorldY));
		F = (1 + GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY));
		G = (1 + GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY));
		H = (1 + GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY));

		float minC = fminf(A, fminf(B, fminf(C, fminf(D, fminf(E, fminf(F, fminf(G, H)))))));
		A = A - (minC - 1);
		B = B - (minC - 1);
		C = C - (minC - 1);
		D = D - (minC - 1);
		E = E - (minC - 1);
		F = F - (minC - 1);
		G = G - (minC - 1);
		H = H - (minC - 1);

		float ABCD, EFGH, ABFE, DCGH, AEHD, BFGC;
		ABCD = (A + B + C + D) / 4;
		EFGH = (E + F + G + H) / 4;
		ABFE = (A + B + F + E) / 4;
		DCGH = (D + C + G + H) / 4;
		AEHD = (A + E + H + D) / 4;
		BFGC = (B + F + G + C) / 4;

		float sum = ABCD + EFGH + ABFE + DCGH + AEHD + BFGC;

		float rand = U.x * sum;
		if (rand < ABCD) {
			Wi = Normalize(Vec3f(U.y, -1, U.z));
			Pdf = ABCD / sum;
		}
		else if (rand < ABCD + EFGH) {
			Wi = Normalize(Vec3f(U.y, 1, U.z));
			Pdf = EFGH / sum;
		}
		else if (rand < ABCD + EFGH + ABFE) {
			Wi = Normalize(Vec3f(-1, U.y, U.z));
			Pdf = ABFE / sum;
		}
		else if (rand < ABCD + EFGH + ABFE + DCGH) {
			Wi = Normalize(Vec3f(1, U.y, U.z));
			Pdf = DCGH / sum;
		}
		else if (rand < ABCD + EFGH + ABFE + DCGH + AEHD) {
			Wi = Normalize(Vec3f(U.y, U.z, -1));
			Pdf = AEHD / sum;
		}
		else {
			Wi = Normalize(Vec3f(U.y, U.z, 1));
			Pdf = BFGC / sum;
		}

		return this->F(Wo, Wi);*/

		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradX, gradY, gradZ;
		gradX = ((D + C + G + H) - (A + B + F + E)) / 4.f;
		gradY = ((E + F + G + H) - (A + B + C + D)) / 4.f;
		gradZ = ((B + F + G + C) - (A + E + H + D)) / 4.f;

		gradX = gradX / 101.f * .99;
		gradY = gradY / 101.f * .99;
		gradZ = gradZ / 101.f * .99;

		float ABCD, EFGH, ABFE, DCGH, AEHD, BFGC;
		ABFE = (1 + gradX);// *(1.f / 6);
		DCGH = (1 - gradX);// * (1.f / 6);
		ABCD = (1 + gradY);// * (1.f / 6);
		EFGH = (1 - gradY);// * (1.f / 6);
		AEHD = (1 + gradZ);// * (1.f / 6);
		BFGC = (1 - gradZ);// * (1.f / 6);

		float rand = U.x * 6;

		float face_probability = -1;
		if (rand < ABCD) {
			Wi = Normalize(Vec3f(U.y, -1, U.z));
			face_probability = ABCD;
			//Pdf = (1.f / 6) / (ABCD * (FOUR_PI_F / 6));
		}
		else if (rand < ABCD + EFGH) {
			Wi = Normalize(Vec3f(U.y, 1, U.z));
			//Pdf = EFGH;
			face_probability = EFGH;
		}
		else if (rand < ABCD + EFGH + ABFE) {
			Wi = Normalize(Vec3f(-1, U.y, U.z));
			//Pdf = ABFE;
			face_probability = ABFE;
		}
		else if (rand < ABCD + EFGH + ABFE + DCGH) {
			Wi = Normalize(Vec3f(1, U.y, U.z));
			//Pdf = DCGH;
			face_probability = DCGH;
		}
		else if (rand < ABCD + EFGH + ABFE + DCGH + AEHD) {
			Wi = Normalize(Vec3f(U.y, U.z, -1));
			//Pdf = AEHD;
			face_probability = AEHD;
		}
		else {
			Wi = Normalize(Vec3f(U.y, U.z, 1));
			//Pdf = BFGC;
			face_probability = BFGC;
		}
		//Pdf = face_probability * INV_4_PI_F;
		Pdf = face_probability * INV_4_PI_F;

		return this->F(Wo, Wi);
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradX, gradY, gradZ;
		gradX = ((D + C + G + H) - (A + B + F + E)) / 4.f;
		gradY = ((E + F + G + H) - (A + B + C + D)) / 4.f;
		gradZ = ((B + F + G + C) - (A + E + H + D)) / 4.f;

		gradX = gradX / 101.f * .99;
		gradY = gradY / 101.f * .99;
		gradZ = gradZ / 101.f * .99;

		float ABCD, EFGH, ABFE, DCGH, AEHD, BFGC;
		ABFE = (1 + gradX);// *(1.f / 6);
		DCGH = (1 - gradX);// * (1.f / 6);
		ABCD = (1 + gradY);// * (1.f / 6);
		EFGH = (1 - gradY);// * (1.f / 6);
		AEHD = (1 + gradZ);// * (1.f / 6);
		BFGC = (1 - gradZ);// * (1.f / 6);

		Vec3f WiAbs(fabsf(Wi.x), fabsf(Wi.y), fabsf(Wi.z));
		if (WiAbs.x == WiAbs.Max()) {
			if (Wi.x > 0)
				return DCGH * INV_4_PI_F;
			else
				return ABFE * INV_4_PI_F;
		}
		else if (WiAbs.y == WiAbs.Max()) {
			if (Wi.y > 0)
				return EFGH * INV_4_PI_F;
			else
				return ABCD * INV_4_PI_F;
		}
		else if (WiAbs.z == WiAbs.Max()) {
			if (Wi.z > 0)
				return BFGC * INV_4_PI_F;
			else
				return AEHD * INV_4_PI_F;
		}

		return 0;
	}

	Vec3f		m_Pe;
	CColorXyz	m_Kd;
};

class CLightPathsOcto
{
public:
	DEV CLightPathsOcto(const Vec3f& Pe, const CColorXyz& Kd) :
		m_Kd(Kd),
		m_Pe(Pe)
	{
	}

	DEV ~CLightPathsOcto(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec3f& U)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		/*float gradX, gradY, gradZ;
		gradX = ((D + C + G + H) / 4.f) - ((A + B + F + E) / 4.f);
		gradY = ((E + F + G + H) / 4.f) - ((A + B + C + D) / 4.f);
		gradZ = ((B + F + G + C) / 4.f) - ((A + E + H + D) / 4.f);

		gradX = gradX / 101.f * .99;
		gradY = gradY / 101.f * .99;
		gradZ = gradZ / 101.f * .99;

		float ABCD, EFGH, ABFE, DCGH, AEHD, BFGC;
		ABFE = (1 + gradX);
		DCGH = (1 - gradX);
		ABCD = (1 + gradY);
		EFGH = (1 - gradY);
		AEHD = (1 + gradZ);
		BFGC = (1 - gradZ);

		float TXnegZneg, TXposZneg, TXposZpos, TXnegZpos, BXnegZneg, BXposZneg, BXposZpos, BXnegZpos;
		TXnegZneg = (EFGH + ABFE + AEHD) / 3.f;
		TXposZneg = (EFGH + DCGH + AEHD) / 3.f;
		TXposZpos = (EFGH + DCGH + BFGC) / 3.f;
		TXnegZpos = (EFGH + ABFE + BFGC) / 3.f;
		BXnegZneg = (ABCD + ABFE + AEHD) / 3.f;
		BXposZneg = (ABCD + DCGH + AEHD) / 3.f;
		BXposZpos = (ABCD + DCGH + BFGC) / 3.f;
		BXnegZpos = (ABCD + ABFE + BFGC) / 3.f;*/

		//float sum = sqrtf(A * A + B * B + C * C + D * D + E * E + F * F + G * G + H * H);
		
		//float min = fminf(A, fminf(B, fminf(C, fminf(D, fminf(E, fminf(F, fminf(G, H))))))) - 1;
		float max = fmaxf(A, fmaxf(B, fmaxf(C, fmaxf(D, fmaxf(E, fmaxf(F, fmaxf(G, H))))))) + 1;
		
		// invert values
		A = max - A;
		B = max - B;
		C = max - C;
		D = max - D;
		E = max - E;
		F = max - F;
		G = max - G;
		H = max - H;

		float sum = A + B + C + D + E + F + G + H;

		float rand = U.x * sum;
		
		float face_probability = -1;

		float S = U.y;
		float T = U.z;
		if (S + T > 1) {
			S = 1 - S;
			T = 1 - T;
		}

		if (rand < E) {
			Vec3f v0 = Vec3f(0, 0, -1);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(-1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = E;
		}
		else if (rand < E + H) {
			Vec3f v0 = Vec3f(1, 0, 0);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, -1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = H;
		}
		else if (rand < E + H + G) {
			Vec3f v0 = Vec3f(0, 0, 1);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = G;
		}
		else if (rand < E + H + G + F) {
			Vec3f v0 = Vec3f(-1, 0, 0);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, 1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = F;
		}
		else if (rand < E + H + G + F + A) {
			Vec3f v0 = Vec3f(0, 0, -1);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(-1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = A;
		}
		else if (rand < E + H + G + F + A + D) {
			Vec3f v0 = Vec3f(1, 0, 0);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, -1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = D;
		}
		else if (rand < E + H + G + F + A + D + C) {
			Vec3f v0 = Vec3f(0, 0, 1);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = C;
		}
		else {
			Vec3f v0 = Vec3f(-1, 0, 0);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, 1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = B;
		}
		//Pdf = (face_probability / length) * INV_4_PI_F;
		//Pdf = (face_probability / length) / (PI_F / 2);
		//Pdf = face_probability * sum (= chance of one face) / area of face (=4*PI/8) -> reduces to (face_probability * 2) / (PI_F * sum)
		Pdf = (face_probability * 2) / (PI_F * sum);

		return this->F(Wo, Wi);
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		/*float gradX, gradY, gradZ;
		gradX = ((D + C + G + H) / 4.f) - ((A + B + F + E) / 4.f);
		gradY = ((E + F + G + H) / 4.f) - ((A + B + C + D) / 4.f);
		gradZ = ((B + F + G + C) / 4.f) - ((A + E + H + D) / 4.f);

		gradX = gradX / 101.f * .99;
		gradY = gradY / 101.f * .99;
		gradZ = gradZ / 101.f * .99;

		float ABCD, EFGH, ABFE, DCGH, AEHD, BFGC;
		ABFE = (1 + gradX);
		DCGH = (1 - gradX);
		ABCD = (1 + gradY);
		EFGH = (1 - gradY);
		AEHD = (1 + gradZ);
		BFGC = (1 - gradZ);

		float TXnegZneg, TXposZneg, TXposZpos, TXnegZpos, BXnegZneg, BXposZneg, BXposZpos, BXnegZpos;
		TXnegZneg = (EFGH + ABFE + AEHD) / 3.f;
		TXposZneg = (EFGH + DCGH + AEHD) / 3.f;
		TXposZpos = (EFGH + DCGH + BFGC) / 3.f;
		TXnegZpos = (EFGH + ABFE + BFGC) / 3.f;
		BXnegZneg = (ABCD + ABFE + AEHD) / 3.f;
		BXposZneg = (ABCD + DCGH + AEHD) / 3.f;
		BXposZpos = (ABCD + DCGH + BFGC) / 3.f;
		BXnegZpos = (ABCD + ABFE + BFGC) / 3.f;
		*/
		
		float max = fmaxf(A, fmaxf(B, fmaxf(C, fmaxf(D, fmaxf(E, fmaxf(F, fmaxf(G, H))))))) + 1;
		float sum = A + B + C + D + E + F + G + H;

		float theta = acosf(Wi.y);
		if (theta >= 0) {
			if (Wi.x <= 0 && Wi.z <= 0)
				return ((max - E) * 2) / (PI_F * sum);
			else if (Wi.x > 0 && Wi.z <= 0)
				return ((max - H) * 2) / (PI_F * sum);
			else if (Wi.x > 0 && Wi.z > 0)
				return ((max - G) * 2) / (PI_F * sum);
			else
				return ((max - F) * 2) / (PI_F * sum);
		}
		else {
			if (Wi.x <= 0 && Wi.z <= 0)
				return ((max - A) * 2) / (PI_F * sum);
			else if (Wi.x > 0 && Wi.z <= 0)
				return ((max - D) * 2) / (PI_F * sum);
			else if (Wi.x > 0 && Wi.z > 0)
				return ((max - C) * 2) / (PI_F * sum);
			else
				return ((max - B) * 2) / (PI_F * sum);
		}
		/*if (theta >= 0) {
			if (Wi.x <= 0 && Wi.z <= 0)
				return ((max - E) * 2) / (PI_F * sum);
			else if (Wi.x > 0 && Wi.z <= 0)
				return (H / length) * INV_4_PI_F;
			else if (Wi.x > 0 && Wi.z > 0)
				return (G / length) * INV_4_PI_F;
			else
				return (F / length) * INV_4_PI_F;
		}
		else {
			if (Wi.x <= 0 && Wi.z <= 0)
				return (A / length) * INV_4_PI_F;
			else if (Wi.x > 0 && Wi.z <= 0)
				return (D / length) * INV_4_PI_F;
			else if (Wi.x > 0 && Wi.z > 0)
				return (C / length) * INV_4_PI_F;
			else
				return (B / length) * INV_4_PI_F;
		}*/

		return 0;
	}

	Vec3f		m_Pe;
	CColorXyz	m_Kd;
};

class CLightPathsOctoGradient
{
public:
	DEV CLightPathsOctoGradient(const Vec3f& Pe, const CColorXyz& Kd) :
		m_Kd(Kd),
		m_Pe(Pe)
	{
	}

	DEV ~CLightPathsOctoGradient(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec3f& U)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = E - C;
		gradBH = H - B;

		gradAG = gradAG / 101.f;// *.99;
		gradDF = gradDF / 101.f;// *.99;
		gradCE = gradCE / 101.f;// *.99;
		gradBH = gradBH / 101.f;// *.99;

		A = (1 + gradAG);
		G = (1 - gradAG);
		D = (1 + gradDF);
		F = (1 - gradDF);
		C = (1 + gradCE);
		E = (1 - gradCE);
		B = (1 + gradBH);
		H = (1 - gradBH);

		float rand = U.x * 8;
		float face_probability = -1;

		float S = U.y;
		float T = U.z;
		if (S + T > 1) {
			S = 1 - S;
			T = 1 - T;
		}

		if (rand < E) {
			Vec3f v0 = Vec3f(0, 0, -1);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(-1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = E;
		}
		else if (rand < E + H) {
			Vec3f v0 = Vec3f(1, 0, 0);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, -1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = H;
		}
		else if (rand < E + H + G) {
			Vec3f v0 = Vec3f(0, 0, 1);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = G;
		}
		else if (rand < E + H + G + F) {
			Vec3f v0 = Vec3f(-1, 0, 0);
			Vec3f v1 = Vec3f(0, 1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, 1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = F;
		}
		else if (rand < E + H + G + F + A) {
			Vec3f v0 = Vec3f(0, 0, -1);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(-1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = A;
		}
		else if (rand < E + H + G + F + A + D) {
			Vec3f v0 = Vec3f(1, 0, 0);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, -1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = D;
		}
		else if (rand < E + H + G + F + A + D + C) {
			Vec3f v0 = Vec3f(0, 0, 1);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(1, 0, 0) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = C;
		}
		else {
			Vec3f v0 = Vec3f(-1, 0, 0);
			Vec3f v1 = Vec3f(0, -1, 0) - v0;
			Vec3f v2 = Vec3f(0, 0, 1) - v0;
			Wi = Normalize(v0 + S * v1 + T * v2);
			face_probability = B;
		}

		//Pdf = (face_probability / 8.0f) * (1.0f / (8 * 4 * PI_F));
		//Pdf = face_probability * 1/8 (=uniform chance of one face) / area of face (=4*PI/8) -> reduces to face_probability * INV_$_PI_F
		//Pdf = (face_probability / 8.0f) / ((4.0f * PI_F) / 8.0f);

		// Since face_probability is between 0 and 2, with 1 represeting uniform chance, this works to adjust throughput
		// to match uniform
		Pdf = face_probability * INV_4_PI_F;

		// Cancel out the inv_4_pi_f on the color F
		//Pdf = INV_4_PI_F;

		// Simply return face_probability so we can use it in LightBalance
		//Pdf = face_probability;

		return this->F(Wo, Wi);
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = C - E;
		gradBH = B - H;

		gradAG = gradAG / 101.f * .99;
		gradDF = gradDF / 101.f * .99;
		gradCE = gradCE / 101.f * .99;
		gradBH = gradBH / 101.f * .99;

		A = (1 + gradAG);
		G = (1 - gradAG);
		D = (1 + gradDF);
		F = (1 - gradDF);
		C = (1 + gradCE);
		E = (1 - gradCE);
		B = (1 + gradBH);
		H = (1 - gradBH);

		float theta = acosf(Wi.y);
		if (theta >= 0) {
			if (Wi.x <= 0 && Wi.z <= 0)
				return E * INV_4_PI_F;
			else if (Wi.x > 0 && Wi.z <= 0)
				return H * INV_4_PI_F;
			else if (Wi.x > 0 && Wi.z > 0)
				return G * INV_4_PI_F;
			else
				return F * INV_4_PI_F;
		}
		else {
			if (Wi.x <= 0 && Wi.z <= 0)
				return A * INV_4_PI_F;
			else if (Wi.x > 0 && Wi.z <= 0)
				return D * INV_4_PI_F;
			else if (Wi.x > 0 && Wi.z > 0)
				return C * INV_4_PI_F;
			else
				return B * INV_4_PI_F;
		}

		return 0;
	}

	Vec3f		m_Pe;
	CColorXyz	m_Kd;
};

class CTestShader
{
public:
	DEV CTestShader(const CColorXyz& Kd) :
		m_Kd(Kd)
	{
	}

	DEV ~CTestShader(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec3f& U)
	{
		Wi = UniformSampleSphere(Vec2f(U.y, U.z));
		//if (U.x <= Chance) {
		//	if (Wi.x >= 0)
		//		Wi.x = Wi.x * -1;
		//}
		//else {
		//	if (Wi.x < 0)
		//		Wi.x = Wi.x * -1;
		//}
		
		Pdf = this->Pdf(Wo, Wi);

		return this->F(Wo, Wi);
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		if (Wi.x <= 0) {
			return Chance / (TWO_PI_F);
		}
		else {
			return (1.f - Chance) / (TWO_PI_F);
		}
	}

	CColorXyz	m_Kd;
	float Chance = 0.5;
};

class CRejectionSampler {
public:
	DEV CRejectionSampler(const Vec3f& Pe, const CColorXyz& Kd) :
		m_Kd(Kd),
		m_Pe(Pe)
	{
	}

	DEV ~CRejectionSampler(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, CBrdfSample& S, CRNG& RNG)
	{
		//return CColorXyz(RNG.Get1());

		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = E - C;
		gradBH = H - B;

		gradAG = gradAG / 102.f;// *.99;
		gradDF = gradDF / 102.f;// *.99;
		gradCE = gradCE / 102.f;// *.99;
		gradBH = gradBH / 102.f;// *.99;

		A = (1.f + gradAG);
		G = (1.f - gradAG);
		D = (1.f + gradDF);
		F = (1.f - gradDF);
		C = (1.f + gradCE);
		E = (1.f - gradCE);
		B = (1.f + gradBH);
		H = (1.f - gradBH);

		//return CColorXyz(G / 2.f, H / 2.f, D / 4.f);

		float AEHD = (A + E + H + D) / 4.f;
		float GCBF = (G + C + B + F) / 4.f;
		float ABFE = (A + B + F + E) / 4.f;
		float GHDC = (G + H + D + C) / 4.f;
		float ABCD = (A + B + C + D) / 4.f;
		float GHEF = (G + H + E + F) / 4.f;
		
		//return CColorXyz(gGradientPower / 2.f);
		
		AEHD = powf(AEHD, gGradientPower);
		GCBF = powf(GCBF, gGradientPower);
		ABFE = powf(ABFE, gGradientPower);
		GHDC = powf(GHDC, gGradientPower);
		ABCD = powf(ABCD, gGradientPower);
		GHEF = powf(GHEF, gGradientPower);
		
		//return CColorXyz(AEHD / 2.f, GCBF / 2.f, ABFE / 2.f);
		//return CColorXyz(GHDC / 2.f, ABCD / 2.f, GHEF / 2.f);
		
		//0.00387597
		/*float faces[8][3] = {
			{GHEF, GCBF, GHDC}, // %G
			{GHEF, ABFE, GCBF}, // %F
			{GHEF, AEHD, ABFE}, // %E
			{GHEF, GHDC, AEHD}, // %H
			{ABCD, GCBF, GHDC}, // %C
			{ABCD, ABFE, GCBF}, // %B
			{ABCD, AEHD, ABFE}, // %A
			{ABCD, GHDC, AEHD}  // %D
		};*/
		float faces[8][3] = {
			{GCBF, GHDC, GHEF}, // %G
			{GCBF, GHEF, ABFE}, // %F
			{GCBF, ABFE, ABCD}, // %E
			{GCBF, ABCD, GHDC}, // %H
			{AEHD, GHDC, GHEF}, // %C
			{AEHD, GHEF, ABFE}, // %B
			{AEHD, ABFE, ABCD}, // %A
			{AEHD, ABCD, GHDC}  // %D
		};

		float sum = 0;
		for (int i = 0; i < 8; i++)
			sum = sum + ((0.5f)*((PI_F - 2.f)*faces[i][0] + faces[i][1] + faces[i][2]));
		CColorXyz Col = CColorXyz(0.5f);
		bool accepted = false;
		int tries = 0;

		//return CColorXyz(sum / 50.f);

		while (!accepted) {
			tries++;

			float rx = RNG.Get1();// S.m_Dir.x;
			float ry = RNG.Get1();// S.m_Dir.y;

			float z = 1.f - 2.f * rx;
			float r = sqrtf(fmaxf(0, 1.f - z * z));
			float phi = acosf(z);
			float theta = 2.f * PI_F * ry;
			float x = r * cosf(theta);
			float y = r * sinf(theta);

			float t_A, t_B, t_C, q;
			if (z > 0) {
				if (theta <= HALF_PI_F) {
					t_A = GCBF;
					t_B = GHDC;
					t_C = GHEF;
					Col = CColorXyz(1, 1, 1);
					q = 0;
				}
				else if (theta <= PI_F) {
					t_A = GCBF;
					t_B = GHEF;
					t_C = ABFE;
					Col = CColorXyz(0, 1, 1);
					q = 1;
				}
				else if (theta <= 3 * HALF_PI_F) {
					t_A = GCBF;
					t_B = ABFE;
					t_C = ABCD;
					Col = CColorXyz(0, 0, 1);
					q = 2;
				}
				else {
					t_A = GCBF;
					t_B = ABCD;
					t_C = GHDC;
					Col = CColorXyz(1, 0, 1);
					q = 3;
				}
			}
			else {
				if (theta <= HALF_PI_F) {
					t_A = AEHD;
					t_B = GHDC;
					t_C = GHEF;
					Col = CColorXyz(1, 1, 0);
					q = 0;
				}
				else if (theta <= PI_F) {
					t_A = AEHD;
					t_B = GHEF;
					t_C = ABFE;
					Col = CColorXyz(0, 1, 0);
					q = 1;
				}
				else if (theta <= 3 * HALF_PI_F) {
					t_A = AEHD;
					t_B = ABFE;
					t_C = ABCD;
					Col = CColorXyz(0, 0, 0);
					q = 2;
				}
				else {
					t_A = AEHD;
					t_B = ABCD;
					t_C = GHDC;
					Col = CColorXyz(1, 0, 0);
					q = 3;
				}
			}

			//float party = fmodf(ry, 0.25) * 4;
			float party = (ry - q * 0.25f) * 4.f;
			//float partx = 1 - fabsf(rx - (1 / 2)) * 2;
			float partx = 1.f - fabsf((phi / PI_F) - 0.5f) * 2.f;
			float hor = (t_C - t_B)*party + t_B;
			float ver = (hor - t_A)*partx + t_A;
			Pdf = ver / sum;
			
			/*
			Vec3f vec = Vec3f(x, y, z);
			float aehd = clamp(Vec3f(0, 0, -1).Dot(vec), 0.f, 1.f) * AEHD;
			float gcbf = clamp(Vec3f(0, 0, 1).Dot(vec), 0.f, 1.f) * GCBF;
			float ghef = clamp(Vec3f(0, 1, 0).Dot(vec), 0.f, 1.f) * GHEF;
			float abcd = clamp(Vec3f(0, -1, 0).Dot(vec), 0.f, 1.f) * ABCD;
			float abfe = clamp(Vec3f(-1, 0, 0).Dot(vec), 0.f, 1.f) * ABFE;
			float ghdc = clamp(Vec3f(1, 0, 0).Dot(vec), 0.f, 1.f) * GHDC;

			float p = aehd + gcbf + ghef + abcd + abfe + ghdc;
			//p /= integral;
			auto a = p / (fabsf(x)+fabsf(y)+fabsf(z)) / sum;

			Pdf = INV_4_PI_F;// / p;
			*/
			//S.LargeStep(RNG);
			//accepted = RNG.Get1() <= Pdf;
			float maxP = fmaxf(AEHD, fmaxf(GCBF, fmaxf(ABFE, fmaxf(GHDC, fmaxf(ABCD, GHEF)))));
			//float c = (maxP / sum) / INV_4_PI_F;
			//accepted = RNG.Get1() <= Pdf / (c*INV_4_PI_F);
			// Due to uniform sample as candidate (https://bookdown.org/rdpeng/advstatcomp/rejection-sampling.html)
			// simplifies to:
			accepted = RNG.Get1() <= ver / maxP;

			//accepted = true;
			//return CColorXyz(t_A, partx, ver / maxP);
			//return CColorXyz((partx + 1.f) / 2.f);
			//return CColorXyz((party + 1.f) / 2.f);
			//return CColorXyz(party);
			//return CColorXyz(phi / PI_F);
			//return CColorXyz((z + 1.f) / 2.f);
			if (accepted) {
				Wi = Vec3f(x, y, z);
				//Wi.Normalize();
				//Pdf = tries;
				
				//if (Wi.y > 0)
				//	return CColorXyz(0);
				//else
				//	return CColorXyz(1.f);
			}
			else {
				//S.LargeStep(RNG);
				/*if (tries > 40) {
					accepted = true;
					Wi = Vec3f(x, y, z);
					Pdf = -1;
				}*/
			}
		}

		//Pdf = this->Pdf(Wo, Wi);

		return this->F(Wo, Wi);
		
		//return CColorXyz(0.25, 0.5, 0.75);
		//return Col;
		/*
		if (Wi.z > 0)
			return CColorXyz(1.f, 0, 0);
		else if (Wi.z < 0)
			return CColorXyz(0.f, 1, 0);
		else
			return CColorXyz(0, 0, 1);
			*/
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = E - C;
		gradBH = H - B;

		gradAG = gradAG / 102.f;// *.99;
		gradDF = gradDF / 102.f;// *.99;
		gradCE = gradCE / 102.f;// *.99;
		gradBH = gradBH / 102.f;// *.99;

		A = (1.f + gradAG);
		G = (1.f - gradAG);
		D = (1.f + gradDF);
		F = (1.f - gradDF);
		C = (1.f + gradCE);
		E = (1.f - gradCE);
		B = (1.f + gradBH);
		H = (1.f - gradBH);

		float AEHD = (A + E + H + D) / 4.f;
		float GCBF = (G + C + B + F) / 4.f;
		float ABFE = (A + B + F + E) / 4.f;
		float GHDC = (G + H + D + C) / 4.f;
		float ABCD = (A + B + C + D) / 4.f;
		float GHEF = (G + H + E + F) / 4.f;

		AEHD = powf(AEHD, gGradientPower);
		GCBF = powf(GCBF, gGradientPower);
		ABFE = powf(ABFE, gGradientPower);
		GHDC = powf(GHDC, gGradientPower);
		ABCD = powf(ABCD, gGradientPower);
		GHEF = powf(GHEF, gGradientPower);

		float faces[8][3] = {
			{GCBF, GHDC, GHEF}, // %G
			{GCBF, GHEF, ABFE}, // %F
			{GCBF, ABFE, ABCD}, // %B
			{GCBF, ABCD, GHDC}, // %C
			{AEHD, GHDC, GHEF}, // %H
			{AEHD, GHEF, ABFE}, // %E
			{AEHD, ABFE, ABCD}, // %A
			{AEHD, ABCD, GHDC}  // %D
		};

		float sum = 0;
		for (int i = 0; i < 8; i++)
			sum = sum + (0.5)*((PI_F - 2.f)*faces[i][0] + faces[i][1] + faces[i][2]);

		float theta = atanf(Wi.y / Wi.x);
		if (Wi.x < 0)
			theta += PI_F;
		else if (Wi.y < 0)
			theta += 2 * PI_F;

		float phi = acosf(Wi.z);

		float t_A, t_B, t_C;
		int face;

		if (Wi.z > 0) {
			if (theta <= HALF_PI_F) {
				face = 0;
			}
			else if (theta <= PI_F) {
				face = 1;
			}
			else if (theta <= 3 * HALF_PI_F) {
				face = 2;
			}
			else {
				face = 3;
			}
		}
		else {
			if (theta <= HALF_PI_F) {
				face = 4;
			}
			else if (theta <= PI_F) {
				face = 5;
			}
			else if (theta <= 3 * HALF_PI_F) {
				face = 6;
			}
			else {
				face = 7;
			}
		}

		t_A = faces[face][2];
		t_B = faces[face][1];
		t_C = faces[face][0];

		float rx = theta / TWO_PI_F;
		float ry = phi / PI_F;

		float partx = (rx - (face % 4) * 0.25f) * 4.f;
		float party = 1.f - fabsf(ry - 0.5f) * 2.f;
		float hor = (t_A - t_B)*partx + t_B;
		float ver = (hor - t_C)*party + t_C;
		float Pdf = ver / sum;

		return Pdf;
	}

	Vec3f		m_Pe;
	CColorXyz	m_Kd;
};

class COneDirectional
{
public:
	HOD COneDirectional(const CColorXyz& Kd) :
		m_Kd(Kd)
	{
	}

	HOD ~COneDirectional(void)
	{
	}

	HOD CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	HOD CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, CRNG& RNG)
	{
		Wi = UniformSampleSphere(RNG.Get2());

		if (RNG.Get1() <= m_p) {
			Wi.z = fabsf(Wi.z);
		}
		else {
			Wi.z = fabsf(Wi.z) * -1.f;
		}
		Pdf = this->Pdf(Wo, Wi);
		return F(Wo, Wi);
	}

	HOD float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		if (Wi.z >= 0) {
			return m_p / TWO_PI_F;
		}
		else {
			return (1.f - m_p) / TWO_PI_F;
		}
	}

	CColorXyz	m_Kd;
	float m_p = 0.25f;
};

class CIsotropicPhase
{
public:
	HOD CIsotropicPhase(const CColorXyz& Kd) :
		m_Kd(Kd)
	{
	}

	HOD ~CIsotropicPhase(void)
	{
	}

	HOD CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	HOD CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec2f& U)
	{
		Wi	= UniformSampleSphere(U);
		Pdf	= this->Pdf(Wo, Wi);

		return F(Wo, Wi);
	}

	HOD float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		return INV_4_PI_F;
	}

	CColorXyz	m_Kd;
};

class CBRDF
{
public:
	HOD CBRDF(const Vec3f& N, const Vec3f& Wo, const CColorXyz& Kd, const CColorXyz& Ks, const float& Ior, const float& Exponent) :
		m_Lambertian(Kd),
		m_Microfacet(Ks, Ior, Exponent),
		m_Nn(N),
		m_Nu(Normalize(Cross(N, Wo))),
		m_Nv(Normalize(Cross(N, m_Nu)))
	{
	}

	HOD ~CBRDF(void)
	{
	}

	HOD Vec3f WorldToLocal(const Vec3f& W)
	{
		return Vec3f(Dot(W, m_Nu), Dot(W, m_Nv), Dot(W, m_Nn));
	}

	HOD Vec3f LocalToWorld(const Vec3f& W)
	{
		return Vec3f(	m_Nu.x * W.x + m_Nv.x * W.y + m_Nn.x * W.z,
						m_Nu.y * W.x + m_Nv.y * W.y + m_Nn.y * W.z,
						m_Nu.z * W.x + m_Nv.z * W.y + m_Nn.z * W.z);
	}

	HOD CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		const Vec3f Wol = WorldToLocal(Wo);
		const Vec3f Wil = WorldToLocal(Wi);

		CColorXyz R;

		R += m_Lambertian.F(Wol, Wil);
		R += m_Microfacet.F(Wol, Wil);

		return R;
	}

	HOD CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const CBrdfSample& S)
	{
		const Vec3f Wol = WorldToLocal(Wo);
		Vec3f Wil;

		CColorXyz R;

		if (S.m_Component <= 0.5f)
		{
			m_Lambertian.SampleF(Wol, Wil, Pdf, S.m_Dir);
		}
		else
		{
			m_Microfacet.SampleF(Wol, Wil, Pdf, S.m_Dir);
		}

		Pdf += m_Lambertian.Pdf(Wol, Wil);
		Pdf += m_Microfacet.Pdf(Wol, Wil);

		R += m_Lambertian.F(Wol, Wil);
		R += m_Microfacet.F(Wol, Wil);

		Wi = LocalToWorld(Wil);

		return R;
	}

	HOD float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		const Vec3f Wol = WorldToLocal(Wo);
		const Vec3f Wil = WorldToLocal(Wi);

		float Pdf = 0.0f;

		Pdf += m_Lambertian.Pdf(Wol, Wil);
		Pdf += m_Microfacet.Pdf(Wol, Wil);

		return Pdf;
	}

	Vec3f			m_Nn;
	Vec3f			m_Nu;
	Vec3f			m_Nv;
	CLambertian		m_Lambertian;
	CMicrofacet		m_Microfacet;
};

class CRejectionSamplerAdvancedFloodfill {
public:
	DEV CRejectionSamplerAdvancedFloodfill(const Vec3f& Pe, const CColorXyz& Kd) :
		m_Kd(Kd),
		m_Pe(Pe)
	{
	}

	DEV ~CRejectionSamplerAdvancedFloodfill(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, CBrdfSample& S, CRNG& RNG)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = E - C;
		gradBH = H - B;

		float maxGradient = fabsf(fmaxf(gradAG, fmaxf(gradDF, fmaxf(gradCE, gradBH))));
		maxGradient = maxGradient == 0.f ? 1.f : maxGradient;
		gradAG = gradAG / maxGradient;
		gradDF = gradDF / maxGradient;
		gradCE = gradCE / maxGradient;
		gradBH = gradBH / maxGradient;

		A = (1.f + gradAG);
		G = (1.f - gradAG);
		D = (1.f + gradDF);
		F = (1.f - gradDF);
		C = (1.f + gradCE);
		E = (1.f - gradCE);
		B = (1.f + gradBH);
		H = (1.f - gradBH);

		float AEHD = (A + E + H + D) / 4.f;
		float GCBF = (G + C + B + F) / 4.f;
		float ABFE = (A + B + F + E) / 4.f;
		float GHDC = (G + H + D + C) / 4.f;
		float ABCD = (A + B + C + D) / 4.f;
		float GHEF = (G + H + E + F) / 4.f;

		AEHD = powf(AEHD, gGradientPower);
		GCBF = powf(GCBF, gGradientPower);
		ABFE = powf(ABFE, gGradientPower);
		GHDC = powf(GHDC, gGradientPower);
		ABCD = powf(ABCD, gGradientPower);
		GHEF = powf(GHEF, gGradientPower);

		float faces[8][3] = {
			{GCBF, GHDC, GHEF}, // %G
			{GCBF, GHEF, ABFE}, // %F
			{GCBF, ABFE, ABCD}, // %E
			{GCBF, ABCD, GHDC}, // %H
			{AEHD, GHDC, GHEF}, // %C
			{AEHD, GHEF, ABFE}, // %B
			{AEHD, ABFE, ABCD}, // %A
			{AEHD, ABCD, GHDC}  // %D
		};

		float sum = 0;
		for (int i = 0; i < 8; i++)
			sum = sum + ((0.5f)*((PI_F - 2.f)*faces[i][0] + faces[i][1] + faces[i][2]));

		bool accepted = false;
		float t_A, t_B, t_C, q, rx, ry, z, r, phi, theta, x, y, party, partx, hor, ver;
		float maxP = fmaxf(AEHD, fmaxf(GCBF, fmaxf(ABFE, fmaxf(GHDC, fmaxf(ABCD, GHEF)))));

		while (!accepted) {
			rx = RNG.Get1();
			ry = RNG.Get1();

			z = 1.f - 2.f * rx;
			r = sqrtf(fmaxf(0, 1.f - z * z));
			phi = acosf(z);
			theta = 2.f * PI_F * ry;
			x = r * cosf(theta);
			y = r * sinf(theta);
			
			if (z > 0) {
				if (theta <= HALF_PI_F) {
					t_A = GCBF;
					t_B = GHDC;
					t_C = GHEF;
					q = 0;
				}
				else if (theta <= PI_F) {
					t_A = GCBF;
					t_B = GHEF;
					t_C = ABFE;
					q = 1;
				}
				else if (theta <= 3 * HALF_PI_F) {
					t_A = GCBF;
					t_B = ABFE;
					t_C = ABCD;
					q = 2;
				}
				else {
					t_A = GCBF;
					t_B = ABCD;
					t_C = GHDC;
					q = 3;
				}
			}
			else {
				if (theta <= HALF_PI_F) {
					t_A = AEHD;
					t_B = GHDC;
					t_C = GHEF;
					q = 0;
				}
				else if (theta <= PI_F) {
					t_A = AEHD;
					t_B = GHEF;
					t_C = ABFE;
					q = 1;
				}
				else if (theta <= 3 * HALF_PI_F) {
					t_A = AEHD;
					t_B = ABFE;
					t_C = ABCD;
					q = 2;
				}
				else {
					t_A = AEHD;
					t_B = ABCD;
					t_C = GHDC;
					q = 3;
				}
			}

			party = (ry - q * 0.25f) * 4.f;
			partx = 1.f - fabsf((phi / PI_F) - 0.5f) * 2.f;
			hor = (t_C - t_B)*party + t_B;
			ver = (hor - t_A)*partx + t_A;
			Pdf = ver / sum;
			
			//float c = (maxP / sum) / INV_4_PI_F;
			//accepted = RNG.Get1() <= Pdf / (c*INV_4_PI_F);
			// Due to uniform sample as candidate (https://bookdown.org/rdpeng/advstatcomp/rejection-sampling.html)
			// simplifies to:
			accepted = RNG.Get1() <= ver / maxP;
			if (accepted) {
				Wi = Vec3f(x, y, z);
			}
		}

		return this->F(Wo, Wi);
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = E - C;
		gradBH = H - B;

		float maxGradient = fabsf(fmaxf(gradAG, fmaxf(gradDF, fmaxf(gradCE, gradBH))));
		maxGradient = maxGradient == 0.f ? 1.f : maxGradient;
		gradAG = gradAG / maxGradient;
		gradDF = gradDF / maxGradient;
		gradCE = gradCE / maxGradient;
		gradBH = gradBH / maxGradient;

		A = (1 + gradAG);
		G = (1 - gradAG);
		D = (1 + gradDF);
		F = (1 - gradDF);
		C = (1 + gradCE);
		E = (1 - gradCE);
		B = (1 + gradBH);
		H = (1 - gradBH);

		float AEHD = (A + E + H + D) / 4;
		float GCBF = (G + C + B + F) / 4;
		float ABFE = (A + B + F + E) / 4;
		float GHDC = (G + H + D + C) / 4;
		float ABCD = (A + B + C + D) / 4;
		float GHEF = (G + H + E + F) / 4;
		
		AEHD = powf(AEHD, gScatteringHeadstart);
		GCBF = powf(GCBF, gScatteringHeadstart);
		ABFE = powf(ABFE, gScatteringHeadstart);
		GHDC = powf(GHDC, gScatteringHeadstart);
		ABCD = powf(ABCD, gScatteringHeadstart);
		GHEF = powf(GHEF, gScatteringHeadstart);

		float faces[8][3] = {
			{GCBF, GHDC, GHEF}, // %G
			{GCBF, GHEF, ABFE}, // %F
			{GCBF, ABFE, ABCD}, // %E
			{GCBF, ABCD, GHDC}, // %H
			{AEHD, GHDC, GHEF}, // %C
			{AEHD, GHEF, ABFE}, // %B
			{AEHD, ABFE, ABCD}, // %A
			{AEHD, ABCD, GHDC}  // %D
		};

		float sum = 0;
		for (int i = 0; i < 8; i++)
			sum = sum + (0.5)*((PI_F - 2.f)*faces[i][0] + faces[i][1] + faces[i][2]);

		//r = x2 + y2 + z2 = 1
		float theta = atanf(Wi.y / Wi.x);
		float phi = acosf(Wi.z);

		float t_A, t_B, t_C, q;
		if (Wi.z > 0) {
			if (theta <= HALF_PI_F) {
				t_A = GCBF;
				t_B = GHDC;
				t_C = GHEF;
				q = 0;
			}
			else if (theta <= PI_F) {
				t_A = GCBF;
				t_B = GHEF;
				t_C = ABFE;
				q = 1;
			}
			else if (theta <= 3 * HALF_PI_F) {
				t_A = GCBF;
				t_B = ABFE;
				t_C = ABCD;
				q = 2;
			}
			else {
				t_A = GCBF;
				t_B = ABCD;
				t_C = GHDC;
				q = 3;
			}
		}
		else {
			if (theta <= HALF_PI_F) {
				t_A = AEHD;
				t_B = GHDC;
				t_C = GHEF;
				q = 0;
			}
			else if (theta <= PI_F) {
				t_A = AEHD;
				t_B = GHEF;
				t_C = ABFE;
				q = 1;
			}
			else if (theta <= 3 * HALF_PI_F) {
				t_A = AEHD;
				t_B = ABFE;
				t_C = ABCD;
				q = 2;
			}
			else {
				t_A = AEHD;
				t_B = ABCD;
				t_C = GHDC;
				q = 3;
			}
		}

		float ry = phi / TWO_PI_F;
		float rx = theta / PI_F;

		float party = (ry - q * 0.25) * 4;
		float partx = 1 - fabsf((phi / PI_F) - 0.5f) * 2;
		float hor = (t_C - t_B)*party + t_B;
		float ver = (hor - t_A)*partx + t_A;
		float Pdf = ver / sum;

		return Pdf;
	}

	Vec3f		m_Pe;
	CColorXyz	m_Kd;
};

class COctoGradientInverseSampler {
public:
	DEV COctoGradientInverseSampler(const Vec3f& Pe, const CColorXyz& Kd) :
		m_Kd(Kd),
		m_Pe(Pe)
	{
	}

	DEV ~COctoGradientInverseSampler(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		return m_Kd * INV_4_PI_F;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const Vec3f& RNG)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = E - C;
		gradBH = H - B;

		gradAG = gradAG / 102.f;// *.99;
		gradDF = gradDF / 102.f;// *.99;
		gradCE = gradCE / 102.f;// *.99;
		gradBH = gradBH / 102.f;// *.99;

		A = (1.f + gradAG);
		G = (1.f - gradAG);
		D = (1.f + gradDF);
		F = (1.f - gradDF);
		C = (1.f + gradCE);
		E = (1.f - gradCE);
		B = (1.f + gradBH);
		H = (1.f - gradBH);

		float AEHD = (A + E + H + D) / 4.f;
		float GCBF = (G + C + B + F) / 4.f;
		float ABFE = (A + B + F + E) / 4.f;
		float GHDC = (G + H + D + C) / 4.f;
		float ABCD = (A + B + C + D) / 4.f;
		float GHEF = (G + H + E + F) / 4.f;

		AEHD = powf(AEHD, gGradientPower);
		GCBF = powf(GCBF, gGradientPower);
		ABFE = powf(ABFE, gGradientPower);
		GHDC = powf(GHDC, gGradientPower);
		ABCD = powf(ABCD, gGradientPower);
		GHEF = powf(GHEF, gGradientPower);

		float faces[8][5] = {
			{GCBF, GHDC, GHEF, 0, 0}, // %G
			{GCBF, GHEF, ABFE, 0, 0}, // %F
			{GCBF, ABFE, ABCD, 0, 0}, // %B
			{GCBF, ABCD, GHDC, 0, 0}, // %C
			{AEHD, GHDC, GHEF, 0, 0}, // %H
			{AEHD, GHEF, ABFE, 0, 0}, // %E
			{AEHD, ABFE, ABCD, 0, 0}, // %A
			{AEHD, ABCD, GHDC, 0, 0}  // %D
		};

		for (int i = 0; i < 8; i++) {
			faces[i][3] = ((0.5f)*((PI_F - 2.f)*faces[i][0] + faces[i][1] + faces[i][2]));
			if (i == 0)
				faces[i][4] = faces[i][3];
			else
				faces[i][4] = faces[i][3] + faces[i - 1][4];
		}

		float rngFace = RNG[0] * faces[7][4];
		float t_A, t_B, t_C, offsetPhi, offsetTheta;
		int face = -1;

		if (rngFace <= faces[0][4]) { // G
			offsetPhi = 1;
			offsetTheta = 0;
			face = 0;
		}
		else if (rngFace <= faces[1][4]) { // F
			offsetPhi = 1;
			offsetTheta = HALF_PI_F;
			face = 1;
		}
		else if (rngFace <= faces[2][4]) { // B
			offsetPhi = 1;
			offsetTheta = PI_F;
			face = 2;
		}
		else if (rngFace <= faces[3][4]) { // C
			offsetPhi = 1;
			offsetTheta = 3.f * HALF_PI_F;
			face = 3;
		}
		else if (rngFace <= faces[4][4]) { // H
			offsetPhi = -1;
			offsetTheta = 0;
			face = 4;
		}
		else if (rngFace <= faces[5][4]) { // E
			offsetPhi = -1;
			offsetTheta = HALF_PI_F;
			face = 5;
		}
		else if (rngFace <= faces[6][4]) { // A
			offsetPhi = -1;
			offsetTheta = PI_F;
			face = 6;
		}
		else if (rngFace <= faces[7][4]) { // D
			offsetPhi = -1;
			offsetTheta = 3.f * HALF_PI_F;
			face = 7;
		}

		t_A = faces[face][2];
		t_B = faces[face][1];
		t_C = faces[face][0];
		float t_D, xInv, yInv, theta, phi, x, y, z, P;

		if (t_A != t_B) {
			float xScaled = RNG[1] * faces[face][3];
			xInv = (PI_F * (sqrtf(8.f * t_A*xScaled + powf(2.f*t_B + (PI_F - 2.f)*t_C, 2.f) - 8.f*t_B*xScaled) - 2.f*t_B - (PI_F - 2.f)*t_C)) / (4.f*(t_A - t_B));
		}
		else {
			xInv = RNG[1] * HALF_PI_F;
		}

		t_D = (t_A - t_B) * (xInv / HALF_PI_F) + t_B;

		if (t_C != t_D) {
			float yScaled = RNG[2] * (t_C + t_D) / 2.f;
			yInv = (t_C - sqrtf(powf(t_C, 2) - 2.f * t_C * yScaled + 2.f * t_D * yScaled)) / (t_C - t_D);
		} else {
			yInv = RNG[2];
		}
		
		theta = xInv + offsetTheta;
		phi = acosf(1.f - yInv);
		x = sinf(phi) * cosf(theta);
		y = sinf(phi) * sinf(theta);
		z = cosf(phi) * offsetPhi;
		Wi = Vec3f(x, y, z);
		
		P = (t_D - t_C)*(phi / HALF_PI_F) + t_C;
		Pdf = P / faces[7][4];

		return this->F(Wo, Wi);
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		float3 mins = floor(make_float3(m_Pe.x, m_Pe.y, m_Pe.z) / gVoxelSizeWorld) * gVoxelSizeWorld;
		float A, B, C, D, E, F, G, H;
		A = GetLightPathValue(mins);
		B = GetLightPathValue(mins + gVoxelSizeWorldZ);
		C = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX);
		D = GetLightPathValue(mins + gVoxelSizeWorldX);

		E = GetLightPathValue(mins + gVoxelSizeWorldY);
		F = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldY);
		G = GetLightPathValue(mins + gVoxelSizeWorldZ + gVoxelSizeWorldX + gVoxelSizeWorldY);
		H = GetLightPathValue(mins + gVoxelSizeWorldX + gVoxelSizeWorldY);

		float gradAG, gradDF, gradCE, gradBH;
		gradAG = G - A;
		gradDF = F - D;
		gradCE = E - C;
		gradBH = H - B;

		gradAG = gradAG / 102.f;// *.99;
		gradDF = gradDF / 102.f;// *.99;
		gradCE = gradCE / 102.f;// *.99;
		gradBH = gradBH / 102.f;// *.99;

		A = (1.f + gradAG);
		G = (1.f - gradAG);
		D = (1.f + gradDF);
		F = (1.f - gradDF);
		C = (1.f + gradCE);
		E = (1.f - gradCE);
		B = (1.f + gradBH);
		H = (1.f - gradBH);

		float AEHD = (A + E + H + D) / 4.f;
		float GCBF = (G + C + B + F) / 4.f;
		float ABFE = (A + B + F + E) / 4.f;
		float GHDC = (G + H + D + C) / 4.f;
		float ABCD = (A + B + C + D) / 4.f;
		float GHEF = (G + H + E + F) / 4.f;

		AEHD = powf(AEHD, gGradientPower);
		GCBF = powf(GCBF, gGradientPower);
		ABFE = powf(ABFE, gGradientPower);
		GHDC = powf(GHDC, gGradientPower);
		ABCD = powf(ABCD, gGradientPower);
		GHEF = powf(GHEF, gGradientPower);

		float faces[8][3] = {
			{GCBF, GHDC, GHEF}, // %G
			{GCBF, GHEF, ABFE}, // %F
			{GCBF, ABFE, ABCD}, // %B
			{GCBF, ABCD, GHDC}, // %C
			{AEHD, GHDC, GHEF}, // %H
			{AEHD, GHEF, ABFE}, // %E
			{AEHD, ABFE, ABCD}, // %A
			{AEHD, ABCD, GHDC}  // %D
		};

		float sum = 0;
		for (int i = 0; i < 8; i++)
			sum = sum + (0.5)*((PI_F - 2.f)*faces[i][0] + faces[i][1] + faces[i][2]);

		float theta = atanf(Wi.y / Wi.x);
		if (Wi.x < 0)
			theta += PI_F;
		else if (Wi.y < 0)
			theta += 2 * PI_F;

		float phi = acosf(Wi.z);

		float t_A, t_B, t_C;
		int face;

		if (Wi.z > 0) {
			if (theta <= HALF_PI_F) {
				face = 0;
			}
			else if (theta <= PI_F) {
				face = 1;
			}
			else if (theta <= 3 * HALF_PI_F) {
				face = 2;
			}
			else {
				face = 3;
			}
		}
		else {
			if (theta <= HALF_PI_F) {
				face = 4;
			}
			else if (theta <= PI_F) {
				face = 5;
			}
			else if (theta <= 3 * HALF_PI_F) {
				face = 6;
			}
			else {
				face = 7;
			}
		}

		t_A = faces[face][2];
		t_B = faces[face][1];
		t_C = faces[face][0];

		float rx = theta / TWO_PI_F;
		float ry = phi / PI_F;

		float partx = (rx - (face % 4) * 0.25f) * 4.f;
		float party = 1.f - fabsf(ry - 0.5f) * 2.f;
		float hor = (t_A - t_B)*partx + t_B;
		float ver = (hor - t_C)*party + t_C;
		float Pdf = ver / sum;

		return Pdf;
	}

	Vec3f		m_Pe;
	CColorXyz	m_Kd;
};

class CVolumeShader
{
public:
	enum EType
	{
		Brdf,
		Phase,
		LightPaths,
		LightPathsOcto,
		LightPathsOctoGradient,
		TestShader,
		LightPathsOctoGradientRejectionSampling,
		OneDirectional,
		LightPathsOctoGradientRejectionSamplingAdvanced,
		OctoGradientInverse

	};

	DEV CVolumeShader(const EType& Type, const Vec3f& Pe, const Vec3f& N, const Vec3f& Wo, const CColorXyz& Kd, const CColorXyz& Ks, const float& Ior, const float& Exponent) :
		m_Type(Type),
		m_Brdf(N, Wo, Kd, Ks, Ior, Exponent),
		m_IsotropicPhase(Kd),
		m_LightPaths(Pe, Kd),
		m_LightPathsOcto(Pe, Kd),
		m_LightPathsOctoGradient(Pe, Kd),
		m_TestShader(Kd),
		m_LightPathsOctoGradientRejectionSampling(Pe, Kd),
		m_OneDirectional(Kd),
		m_RejectionAdvancedFloodfill(Pe, Kd),
		m_OctoGradientInverse(Pe, Kd)
	{
	}

	DEV ~CVolumeShader(void)
	{
	}

	DEV CColorXyz F(const Vec3f& Wo, const Vec3f& Wi)
	{
		switch (m_Type)
		{
			case Brdf:
				return m_Brdf.F(Wo, Wi);

			case Phase:
				return m_IsotropicPhase.F(Wo, Wi);

			case LightPaths:
				return m_LightPaths.F(Wo, Wi);

			case LightPathsOcto:
				return m_LightPathsOcto.F(Wo, Wi);

			case LightPathsOctoGradient:
				return m_LightPathsOctoGradient.F(Wo, Wi);

			case TestShader:
				return m_TestShader.F(Wo, Wi);

			case LightPathsOctoGradientRejectionSampling:
				return m_LightPathsOctoGradientRejectionSampling.F(Wo, Wi);

			case OneDirectional:
				return m_OneDirectional.F(Wo, Wi);

			case LightPathsOctoGradientRejectionSamplingAdvanced:
				return m_RejectionAdvancedFloodfill.F(Wo, Wi);
			case OctoGradientInverse:
				return m_OctoGradientInverse.F(Wo, Wi);
		}

		return 1.0f;
	}

	DEV CColorXyz SampleF(const Vec3f& Wo, Vec3f& Wi, float& Pdf, const CBrdfSample& S)
	{
		switch (m_Type)
		{
			case Brdf:
				return m_Brdf.SampleF(Wo, Wi, Pdf, S);

			case Phase:
				return m_IsotropicPhase.SampleF(Wo, Wi, Pdf, S.m_Dir);

			case LightPaths:
				return m_LightPaths.SampleF(Wo, Wi, Pdf, Vec3f(S.m_Component, S.m_Dir.x, S.m_Dir.y));

			case LightPathsOcto:
				return m_LightPathsOcto.SampleF(Wo, Wi, Pdf, Vec3f(S.m_Component, S.m_Dir.x, S.m_Dir.y));

			case LightPathsOctoGradient:
				return m_LightPathsOctoGradient.SampleF(Wo, Wi, Pdf, Vec3f(S.m_Component, S.m_Dir.x, S.m_Dir.y));

			case TestShader:
				return m_TestShader.SampleF(Wo, Wi, Pdf, Vec3f(S.m_Component, S.m_Dir.x, S.m_Dir.y));

			case OctoGradientInverse:
				return m_OctoGradientInverse.SampleF(Wo, Wi, Pdf, Vec3f(S.m_Component, S.m_Dir.x, S.m_Dir.y));
		}
			
	}

	DEV CColorXyz SampleFRejection(const Vec3f& Wo, Vec3f& Wi, float& Pdf, CBrdfSample& S, CRNG& RNG) {
		switch (m_Type)
		{
			case LightPathsOctoGradientRejectionSampling:
				return m_LightPathsOctoGradientRejectionSampling.SampleF(Wo, Wi, Pdf, S, RNG);

			case OneDirectional:
				return m_OneDirectional.SampleF(Wo, Wi, Pdf, RNG);

			case LightPathsOctoGradientRejectionSamplingAdvanced:
				return m_RejectionAdvancedFloodfill.SampleF(Wo, Wi, Pdf, S, RNG);
		}
	}

	DEV float Pdf(const Vec3f& Wo, const Vec3f& Wi)
	{
		switch (m_Type)
		{
			case Brdf:
				return m_Brdf.Pdf(Wo, Wi);

			case Phase:
				return m_IsotropicPhase.Pdf(Wo, Wi);

			case LightPaths:
				return m_LightPaths.Pdf(Wo, Wi);

			case LightPathsOcto:
				return m_LightPathsOcto.Pdf(Wo, Wi);

			case LightPathsOctoGradient:
				return m_LightPathsOctoGradient.Pdf(Wo, Wi);

			case TestShader:
				return m_TestShader.Pdf(Wo, Wi);

			case LightPathsOctoGradientRejectionSampling:
				return m_LightPathsOctoGradientRejectionSampling.Pdf(Wo, Wi);

			case OneDirectional:
				return m_OneDirectional.Pdf(Wo, Wi);

			case LightPathsOctoGradientRejectionSamplingAdvanced:
				return m_RejectionAdvancedFloodfill.Pdf(Wo, Wi);

			case OctoGradientInverse:
				return m_OctoGradientInverse.Pdf(Wo, Wi);
		}

		return 1.0f;
	}

	EType				m_Type;
	CBRDF				m_Brdf;
	CIsotropicPhase		m_IsotropicPhase;
	CLightPaths			m_LightPaths;
	CLightPathsOcto		m_LightPathsOcto;
	CLightPathsOctoGradient m_LightPathsOctoGradient;
	CTestShader			m_TestShader;
	CRejectionSampler	m_LightPathsOctoGradientRejectionSampling;
	COneDirectional		m_OneDirectional;
	CRejectionSamplerAdvancedFloodfill m_RejectionAdvancedFloodfill;
	COctoGradientInverseSampler m_OctoGradientInverse;
};