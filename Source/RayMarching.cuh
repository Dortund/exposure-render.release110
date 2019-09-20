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
#include "Scene.h"
#include "CudaUtilities.h"
#include <curand_kernel.h>

#define KRNL_SS_BLOCK_W		16
#define KRNL_SS_BLOCK_H		8
#define KRNL_SS_BLOCK_SIZE	KRNL_SS_BLOCK_W * KRNL_SS_BLOCK_H

DEV inline bool SampleDistanceRMpropertyBased(CRay& R, CRNG& RNG, Vec3f& Ps)
{
	const int TID = threadIdx.y * blockDim.x + threadIdx.x;

	__shared__ float MinT[KRNL_SS_BLOCK_SIZE];
	__shared__ float MaxT[KRNL_SS_BLOCK_SIZE];

	// Check if the ray intersects the volume's bounding box. If not, there's no use calculating anything
	if (!IntersectBox(R, &MinT[TID], &MaxT[TID]))
		return false;

	// Upated values to stay within bounds specified in the ray.
	MinT[TID] = max(MinT[TID], R.m_MinT);
	MaxT[TID] = min(MaxT[TID], R.m_MaxT);


	const float S = -log(RNG.Get1()) / gDensityScale;
	//const float S = 1.f;
	float Sum = 0.0f;
	float SigmaT = 0.0f;

	// Increase MinT by a random factor of the stepsize, to create the first POSSIBLE scattering point
	MinT[TID] += RNG.Get1() * gStepSize;
	//int steps = -1;
	//float maxOp = -1;
	while (Sum < S)
	{
		//steps++;
		// Scattering point is MinT (distance) * ray direction from ray origin
		Ps = R.m_O + MinT[TID] * R.m_D;

		// If MinT exceeds MaxT, where out of the volume and should stop looking for more points
		if (MinT[TID] > MaxT[TID])
			return false;
			//return Vec3f(S, Sum, maxOp);

		// Calculate SigmaT
		// TODO What exactly is SigmaT???
		materialProperties properties;
		sampleProperties(properties, nullptr, nullptr, Ps);
		//SigmaT = gDensityScale * properties.opacity;
		SigmaT = properties.opacity;

		//maxOp = Fmaxf(maxOp, properties.opacity);

		// Increase Sum by the calculated SigmaT * Stepsize
		Sum += SigmaT * (gStepSize / fminf(gVoxelSizeWorld.x, fminf(gVoxelSizeWorld.y, gVoxelSizeWorld.z)));
		
		// Increase MinT so it can be used to calculate the next possible scattering point in the next loop
		MinT[TID] += gStepSize;
	}

	//return Vec3f(S, Sum, maxOp);
	return true;
}

DEV inline bool SampleDistanceRM(CRay& R, curandState state, Vec3f& Ps) {
	const int TID = threadIdx.y * blockDim.x + threadIdx.x;

	__shared__ float MinT[KRNL_SS_BLOCK_SIZE];
	__shared__ float MaxT[KRNL_SS_BLOCK_SIZE];

	// Check if the ray intersects the volume's bounding box. If not, there's no use calculating anything
	if (!IntersectBox(R, &MinT[TID], &MaxT[TID]))
		return false;

	// Upated values to stay within bounds specified in the ray.
	MinT[TID] = max(MinT[TID], R.m_MinT);
	MaxT[TID] = min(MaxT[TID], R.m_MaxT);


	const float S = -log(curand_uniform(&state)) / gDensityScale;
	float Sum = 0.0f;
	float SigmaT = 0.0f;

	// Increase MinT by a random factor of the stepsize, to create the first POSSIBLE scattering point
	MinT[TID] += curand_uniform(&state) * gStepSize;

	while (Sum < S)
	{
		// Scattering point is MinT (distance) * ray direction from ray origin
		Ps = R.m_O + MinT[TID] * R.m_D;

		// If MinT exceeds MaxT, where out of the volume and should stop looking for more points
		if (MinT[TID] > MaxT[TID])
			return false;

		// Calculate SigmaT
		// TODO What exactly is SigmaT???
		SigmaT = gDensityScale * GetOpacity(GetNormalizedIntensity(Ps));

		// Increase Sum by the calculated SigmaT * Stepsize
		Sum += SigmaT * gStepSize;
		// Increase MinT so it can be used to calculate the next possible scattering point in the next loop
		MinT[TID] += gStepSize;
	}

	return true;
}

/// <summary>
/// Calculates the position of the scattering point using ???Woodcock Tracking???
/// </summary>
/// <param name="Ps">Pointer to variable to hold the position of the scattering point</param>
DEV inline bool SampleDistanceRM(CRay& R, CRNG& RNG, Vec3f& Ps)
{
	const int TID = threadIdx.y * blockDim.x + threadIdx.x;

	__shared__ float MinT[KRNL_SS_BLOCK_SIZE];
	__shared__ float MaxT[KRNL_SS_BLOCK_SIZE];

	// Check if the ray intersects the volume's bounding box. If not, there's no use calculating anything
	if (!IntersectBox(R, &MinT[TID], &MaxT[TID]))
		return false;

	// Upated values to stay within bounds specified in the ray.
	MinT[TID] = max(MinT[TID], R.m_MinT);
	MaxT[TID] = min(MaxT[TID], R.m_MaxT);


	const float S	= -log(RNG.Get1()) / gDensityScale;
	float Sum		= 0.0f;
	float SigmaT	= 0.0f;

	// Increase MinT by a random factor of the stepsize, to create the first POSSIBLE scattering point
	MinT[TID] += RNG.Get1() * gStepSize;

	while (Sum < S)
	{
		// Scattering point is MinT (distance) * ray direction from ray origin
		Ps = R.m_O + MinT[TID] * R.m_D;

		// If MinT exceeds MaxT, where out of the volume and should stop looking for more points
		if (MinT[TID] > MaxT[TID])
			return false;
		
		// Calculate SigmaT
		// TODO What exactly is SigmaT???
		SigmaT	= gDensityScale * GetOpacity(GetNormalizedIntensity(Ps));

		// Increase Sum by the calculated SigmaT * Stepsize
		Sum			+= SigmaT * gStepSize;
		// Increase MinT so it can be used to calculate the next possible scattering point in the next loop
		MinT[TID]	+= gStepSize;
	}

	return true;
}

DEV inline bool FreePathRMPropertyBased(CRay& R, CRNG& RNG)
{
	const int TID = threadIdx.y * blockDim.x + threadIdx.x;

	__shared__ float MinT[KRNL_SS_BLOCK_SIZE];
	__shared__ float MaxT[KRNL_SS_BLOCK_SIZE];
	__shared__ Vec3f Ps[KRNL_SS_BLOCK_SIZE];

	if (!IntersectBox(R, &MinT[TID], &MaxT[TID]))
		return false;

	MinT[TID] = max(MinT[TID], R.m_MinT);
	MaxT[TID] = min(MaxT[TID], R.m_MaxT);

	const float S = -log(RNG.Get1()) / gDensityScale;
	float Sum = 0.0f;
	float SigmaT = 0.0f;

	MinT[TID] += RNG.Get1() * gStepSizeShadow;

	while (Sum < S)
	{
		Ps[TID] = R.m_O + MinT[TID] * R.m_D;

		if (MinT[TID] > MaxT[TID])
			return false;

		materialProperties properties;
		sampleProperties(properties, nullptr, nullptr, Ps[TID]);
		//SigmaT = gDensityScale * properties.opacity;
		SigmaT = properties.opacity;

		Sum += SigmaT * gStepSizeShadow;
		MinT[TID] += gStepSizeShadow;
	}

	return true;
}

DEV inline bool FreePathRM(CRay& R, CRNG& RNG)
{
	const int TID = threadIdx.y * blockDim.x + threadIdx.x;

	__shared__ float MinT[KRNL_SS_BLOCK_SIZE];
	__shared__ float MaxT[KRNL_SS_BLOCK_SIZE];
	__shared__ Vec3f Ps[KRNL_SS_BLOCK_SIZE];

	if (!IntersectBox(R, &MinT[TID], &MaxT[TID]))
		return false;

	MinT[TID] = max(MinT[TID], R.m_MinT);
	MaxT[TID] = min(MaxT[TID], R.m_MaxT);

	const float S	= -log(RNG.Get1()) / gDensityScale;
	float Sum		= 0.0f;
	float SigmaT	= 0.0f;

	MinT[TID] += RNG.Get1() * gStepSizeShadow;
	
	while (Sum < S)
	{
		Ps[TID] = R.m_O + MinT[TID] * R.m_D;

		if (MinT[TID] > MaxT[TID])
			return false;
		
		SigmaT = gDensityScale * GetOpacity(GetNormalizedIntensity(Ps[TID]));

		Sum			+= SigmaT * gStepSizeShadow;
		MinT[TID]	+= gStepSizeShadow;
	}

	return true;
}

DEV inline bool NearestIntersection(CRay R, CScene* pScene, float& T)
{
	float MinT = 0.0f, MaxT = 0.0f;

	if (!IntersectBox(R, &MinT, &MaxT))
		return false;

	MinT = max(MinT, R.m_MinT);
	MaxT = min(MaxT, R.m_MaxT);

	Vec3f Ps; 

	T = MinT;

	while (T < MaxT)
	{
		Ps = R.m_O + T * R.m_D;

		if (GetOpacity(GetNormalizedIntensity(Ps)) > 0.0f)
			return true;

		T += gStepSize;
	}

	return false;
}