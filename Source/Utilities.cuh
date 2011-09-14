#pragma once

#include "Scene.h"

__device__ inline float Dens(CScene* pScene, const Vec3f& P)
{
	return ((float)(SHRT_MAX) * tex3D(gTexDensity, P.x / pScene->m_BoundingBox.LengthX(), P.y / pScene->m_BoundingBox.LengthY(), P.z /  pScene->m_BoundingBox.LengthZ()));
}

__device__ inline Vec3f NormalizedGradient(CScene* pScene, const Vec3f& P)
{
	Vec3f Normal;

	float Delta = 0.01f;

	Normal.x = Dens(pScene, P + Vec3f(Delta, 0.0f, 0.0f)) - Dens(pScene, P - Vec3f(Delta, 0.0f, 0.0f));
	Normal.y = Dens(pScene, P + Vec3f(0.0f, Delta, 0.0f)) - Dens(pScene, P - Vec3f(0.0f, Delta, 0.0f));
	Normal.z = Dens(pScene, P + Vec3f(0.0f, 0.0f, Delta)) - Dens(pScene, P - Vec3f(0.0f, 0.0f, Delta));

	Normal.Normalize();

	return -Normal;
}

DEV float Density(CScene* pScene, const Vec3f& P)
{
	return ((float)(SHRT_MAX) * tex3D(gTexDensity, P.x / pScene->m_BoundingBox.LengthX(), P.y / pScene->m_BoundingBox.LengthY(), P.z /  pScene->m_BoundingBox.LengthZ()));
}

DEV float Extinction(CScene* pScene, const Vec3f& P)
{
	return tex3D(gTexExtinction, P.x / pScene->m_BoundingBox.LengthX(), P.y / pScene->m_BoundingBox.LengthY(), P.z /  pScene->m_BoundingBox.LengthZ());
}

DEV CColorRgbHdr GetOpacity(CScene* pScene, const float& D)
{
	return pScene->m_DensityScale * pScene->m_TransferFunctions.m_Opacity.F(D);
}

DEV CColorRgbHdr GetDiffuse(CScene* pScene, const float& D)
{
	return pScene->m_TransferFunctions.m_Diffuse.F(D);
}

DEV CColorRgbHdr GetSpecular(CScene* pScene, const float& D)
{
	return pScene->m_TransferFunctions.m_Specular.F(D);
}

DEV CColorRgbHdr GetEmission(CScene* pScene, const float& D)
{
	return pScene->m_TransferFunctions.m_Emission.F(D);
}

DEV CColorRgbHdr GetRoughness(CScene* pScene, const float& D)
{
	return pScene->m_TransferFunctions.m_Roughness.F(D);
}

DEV bool NearestLight(CScene* pScene, CRay& R, CColorXyz& LightColor)
{
	// Whether a hit with a light was found or not 
	bool Hit = false;
	
	float T = 0.0f;

	CRay RayCopy = R;

	for (int i = 0; i < pScene->m_Lighting.m_NoLights; i++)
	{
		if (pScene->m_Lighting.m_Lights[i].Intersect(RayCopy, T, LightColor))
			Hit = true;
	}
	
	return Hit;
}