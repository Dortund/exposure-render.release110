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
#include "Variance.h"

#define KRNL_TM_BLOCK_W		8
#define KRNL_TM_BLOCK_H		8
#define KRNL_TM_BLOCK_SIZE	KRNL_TM_BLOCK_W * KRNL_TM_BLOCK_H

KERNEL void KrnlToneMap(CCudaView* pView)
{
	const int X 	= blockIdx.x * blockDim.x + threadIdx.x;
	const int Y		= blockIdx.y * blockDim.y + threadIdx.y;

	if (X >= gFilmWidth || Y >= gFilmHeight)
		return;

	const CColorXyza Color = pView->m_RunningEstimateXyza.Get(X, Y);

	CColorRgbHdr RgbHdr;

	RgbHdr.FromXYZ(Color.c[0], Color.c[1], Color.c[2]);

	RgbHdr.r = Clamp(1.0f - expf(-(RgbHdr.r * gInvExposure)), 0.0, 1.0f);
	RgbHdr.g = Clamp(1.0f - expf(-(RgbHdr.g * gInvExposure)), 0.0, 1.0f);
	RgbHdr.b = Clamp(1.0f - expf(-(RgbHdr.b * gInvExposure)), 0.0, 1.0f);

	if (Color.c[1] <= 1000)
		pView->m_EstimateRgbaLdr.Set(CColorRgbaLdr(RgbHdr.r * 255.0f, RgbHdr.g * 255.0f, RgbHdr.b * 255.0f, 0), X, Y);
	else
		pView->m_EstimateRgbaLdr.Set(CColorRgbaLdr(255.f, 0, 0, 0), X, Y);
}

void ToneMap(CScene* pScene, CScene* pDevScene, CCudaView* pDevView)
{
	const dim3 KernelBlock(KRNL_TM_BLOCK_W, KRNL_TM_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.GetWidth() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.GetHeight() / (float)KernelBlock.y));

	KrnlToneMap<<<KernelGrid, KernelBlock>>>(pDevView);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Tone Map");
}

KERNEL void KrnlToneMapCopy(CCudaView* pView)
{
	const int X = blockIdx.x * blockDim.x + threadIdx.x;
	const int Y = blockIdx.y * blockDim.y + threadIdx.y;

	if (X >= gFilmWidth || Y >= gFilmHeight)
		return;

	const CColorXyza Color = pView->m_RunningEstimateXyza.Get(X, Y);

	CColorRgbHdr RgbHdr = CColorRgbHdr(Color.c[0], Color.c[1], Color.c[2]);

	//RgbHdr.FromXYZ(Color.c[0], Color.c[1], Color.c[2]);


	RgbHdr.r = Clamp(RgbHdr.r * gInvExposure, 0.0, 1.0f);
	RgbHdr.g = Clamp(RgbHdr.g * gInvExposure, 0.0, 1.0f);
	RgbHdr.b = Clamp(RgbHdr.b * gInvExposure, 0.0, 1.0f);
	Vec3f dir = Vec3f(RgbHdr.r, RgbHdr.g, RgbHdr.b);
	dir = dir * 2 - 1;
	dir.Normalize();
	pView->m_EstimateRgbaLdr.Set(CColorRgbaLdr(RgbHdr.r * 255.0f, RgbHdr.g * 255.0f, RgbHdr.b * 255.0f, 0), X, Y);
	/*pView->m_EstimateRgbaLdr.Set(CColorRgbaLdr(
		(dir.x + 1.f) / 2.f * 255.0f,
		(dir.y + 1.f) / 2.f * 255.0f,
		(dir.z + 1.f) / 2.f * 255.0f,
		0), X, Y);
		*/
}

void ToneMapCopy(CScene* pScene, CScene* pDevScene, CCudaView* pDevView)
{
	const dim3 KernelBlock(KRNL_TM_BLOCK_W, KRNL_TM_BLOCK_H);
	const dim3 KernelGrid((int)ceilf((float)pScene->m_Camera.m_Film.GetWidth() / (float)KernelBlock.x), (int)ceilf((float)pScene->m_Camera.m_Film.GetHeight() / (float)KernelBlock.y));

	KrnlToneMapCopy<<<KernelGrid, KernelBlock>>>(pDevView);
	cudaThreadSynchronize();
	HandleCudaKernelError(cudaGetLastError(), "Tone Map Empty");
}