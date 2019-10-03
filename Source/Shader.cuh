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
			Vec3f v0 = Vec3f(0, 0, 1);
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
			Vec3f v0 = Vec3f(0, 0, 1);
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

class CVolumeShader
{
public:
	enum EType
	{
		Brdf,
		Phase,
		LightPaths,
		LightPathsOcto
	};

	DEV CVolumeShader(const EType& Type, const Vec3f& Pe, const Vec3f& N, const Vec3f& Wo, const CColorXyz& Kd, const CColorXyz& Ks, const float& Ior, const float& Exponent) :
		m_Type(Type),
		m_Brdf(N, Wo, Kd, Ks, Ior, Exponent),
		m_IsotropicPhase(Kd),
		m_LightPaths(Pe, Kd),
		m_LightPathsOcto(Pe, Kd)
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
		}

		return 1.0f;
	}

	EType				m_Type;
	CBRDF				m_Brdf;
	CIsotropicPhase		m_IsotropicPhase;
	CLightPaths			m_LightPaths;
	CLightPathsOcto		m_LightPathsOcto;
};