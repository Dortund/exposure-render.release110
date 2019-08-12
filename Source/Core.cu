/*
	Copyright (c) 2011, T. Kroes <t.kroes@tudelft.nl>
	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

	- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
	- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
	- Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
	
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "Core.cuh"

texture<short, cudaTextureType3D, cudaReadModeNormalizedFloat>		gTexDensity;
texture<short, cudaTextureType3D, cudaReadModeNormalizedFloat>		gTexGradientMagnitude;
texture<float, cudaTextureType3D, cudaReadModeElementType>			gTexExtinction;
texture<float, cudaTextureType1D, cudaReadModeElementType>			gTexOpacity;
texture<float4, cudaTextureType1D, cudaReadModeElementType>			gTexDiffuse;
texture<float4, cudaTextureType1D, cudaReadModeElementType>			gTexSpecular;
texture<float, cudaTextureType1D, cudaReadModeElementType>			gTexRoughness;
texture<float4, cudaTextureType1D, cudaReadModeElementType>			gTexEmission;
texture<uchar4, cudaTextureType2D, cudaReadModeNormalizedFloat>		gTexRunningEstimateRgba;
texture<float4, cudaTextureType3D, cudaReadModeElementType>			gTexOpacityGradient;
texture<float, cudaTextureType3D, cudaReadModeNormalizedFloat>		gTexOpacityMagnitudeNormalized;

cudaArray* gpDensityArray				= NULL;
cudaArray* gpGradientMagnitudeArray		= NULL;
cudaArray* gpOpacityArray				= NULL;
cudaArray* gpDiffuseArray				= NULL;
cudaArray* gpSpecularArray				= NULL;
cudaArray* gpRoughnessArray				= NULL;
cudaArray* gpEmissionArray				= NULL;

CD float3		gAaBbMin;
CD float3		gAaBbMax;
CD float3		gInvAaBbMin;
CD float3		gInvAaBbMax;
CD float		gIntensityMin;
CD float		gIntensityMax;
CD float		gIntensityRange;
CD float		gIntensityInvRange;
CD float		gStepSize;
CD float		gStepSizeShadow;
CD float		gDensityScale;
CD float		gGradientDelta;
CD float		gInvGradientDelta;
CD float3		gGradientDeltaX;
CD float3		gGradientDeltaY;
CD float3		gGradientDeltaZ;
CD int			gFilmWidth;
CD int			gFilmHeight;
CD int			gFilmNoPixels;
CD int			gFilterWidth;
CD float		gFilterWeights[10];
CD float		gExposure;
CD float		gInvExposure;
CD float		gGamma;
CD float		gInvGamma;
CD float		gDenoiseEnabled;
CD float		gDenoiseWindowRadius;
CD float		gDenoiseInvWindowArea;
CD float		gDenoiseNoise;
CD float		gDenoiseWeightThreshold;
CD float		gDenoiseLerpThreshold;
CD float		gDenoiseLerpC;
CD float		gNoIterations;
CD float		gInvNoIterations;
CD float		gScatteringHeadstart;

#define TF_NO_SAMPLES		128
#define INV_TF_NO_SAMPLES	1.0f / (float)TF_NO_SAMPLES

#include "Model.cuh"
#include "View.cuh"
#include "Blur.cuh"
#include "Denoise.cuh"
#include "Estimate.cuh"
#include "Utilities.cuh"
#include "SingleScattering.cuh"
#include "MultipleScattering.cuh"
#include "NearestIntersection.cuh"
#include "SpecularBloom.cuh"
#include "ToneMap.cuh"
#include "PreCalculatedScattering.cuh"
#include <iostream>

CCudaModel	gModel;
CCudaView	gRenderCanvasView;
CCudaView	gNavigatorView;

// pre-calculated
int m_nrPoints = -1;
Vec3f* m_DevPoints = NULL;
CColorXyza* m_DevColours = NULL;
float* m_DevConnections = NULL;

void BindDensityBuffer(short* pBuffer, cudaExtent Extent)
{
	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<short>();

	HandleCudaError(cudaMalloc3DArray(&gpDensityArray, &ChannelDesc, Extent));

	cudaMemcpy3DParms CopyParams = {0};

	CopyParams.srcPtr	= make_cudaPitchedPtr(pBuffer, Extent.width * sizeof(short), Extent.width, Extent.height);
	CopyParams.dstArray	= gpDensityArray;
	CopyParams.extent	= Extent;
	CopyParams.kind		= cudaMemcpyHostToDevice;
	
	HandleCudaError(cudaMemcpy3D(&CopyParams));

	gTexDensity.normalized		= true;
	gTexDensity.filterMode		= cudaFilterModeLinear;      
	gTexDensity.addressMode[0]	= cudaAddressModeClamp;  
	gTexDensity.addressMode[1]	= cudaAddressModeClamp;
  	gTexDensity.addressMode[2]	= cudaAddressModeClamp;

	HandleCudaError(cudaBindTextureToArray(gTexDensity, gpDensityArray, ChannelDesc));
}

void BindGradientMagnitudeBuffer(short* pBuffer, cudaExtent Extent)
{
	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<short>();
	HandleCudaError(cudaMalloc3DArray(&gpGradientMagnitudeArray, &ChannelDesc, Extent));

	cudaMemcpy3DParms CopyParams = {0};

	CopyParams.srcPtr	= make_cudaPitchedPtr(pBuffer, Extent.width * sizeof(short), Extent.width, Extent.height);
	CopyParams.dstArray	= gpGradientMagnitudeArray;
	CopyParams.extent	= Extent;
	CopyParams.kind		= cudaMemcpyHostToDevice;
	
	HandleCudaError(cudaMemcpy3D(&CopyParams));

	gTexGradientMagnitude.normalized		= true;
	gTexGradientMagnitude.filterMode		= cudaFilterModeLinear;      
	gTexGradientMagnitude.addressMode[0]	= cudaAddressModeClamp;  
	gTexGradientMagnitude.addressMode[1]	= cudaAddressModeClamp;
  	gTexGradientMagnitude.addressMode[2]	= cudaAddressModeClamp;

	HandleCudaError(cudaBindTextureToArray(gTexGradientMagnitude, gpGradientMagnitudeArray, ChannelDesc));
}

void UnbindDensityBuffer(void)
{
	HandleCudaError(cudaFreeArray(gpDensityArray));
	gpDensityArray = NULL;
	HandleCudaError(cudaUnbindTexture(gTexDensity));
}

void UnbindGradientMagnitudeBuffer(void)
{
	HandleCudaError(cudaFreeArray(gpGradientMagnitudeArray));
	gpGradientMagnitudeArray = NULL;
	HandleCudaError(cudaUnbindTexture(gTexGradientMagnitude));
}

void BindRenderCanvasView(const CResolution2D& Resolution)
{
	gRenderCanvasView.Resize(Resolution);

	cudaChannelFormatDesc Channel;
	
	Channel = cudaCreateChannelDesc<uchar4>();

	HandleCudaError(cudaBindTexture2D(0, gTexRunningEstimateRgba, gRenderCanvasView.m_EstimateRgbaLdr.GetPtr(), Channel, gRenderCanvasView.GetWidth(), gRenderCanvasView.GetHeight(), gRenderCanvasView.m_EstimateRgbaLdr.GetPitch()));
}

void ResetRenderCanvasView(void)
{
	gRenderCanvasView.Reset();
}

void FreeRenderCanvasView(void)
{
	gRenderCanvasView.Free();
}

unsigned char* GetDisplayEstimate(void)
{
	return (unsigned char*)gRenderCanvasView.m_DisplayEstimateRgbLdr.GetPtr(0, 0);
}

unsigned char* GetFrameDisplayEstimate(void)
{
	return (unsigned char*)gRenderCanvasView.m_FrameEstimateXyza.GetPtr(0, 0);
}

void BindTransferFunctionOpacity(CTransferFunction& TransferFunctionOpacity)
{
	gTexOpacity.normalized		= true;
	gTexOpacity.filterMode		= cudaFilterModeLinear;
	gTexOpacity.addressMode[0]	= cudaAddressModeClamp;

	float Opacity[TF_NO_SAMPLES];

	for (int i = 0; i < TF_NO_SAMPLES; i++)
		Opacity[i] = TransferFunctionOpacity.F((float)i * INV_TF_NO_SAMPLES).r;
	
	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<float>();

	if (gpOpacityArray == NULL)
		HandleCudaError(cudaMallocArray(&gpOpacityArray, &ChannelDesc, TF_NO_SAMPLES, 1));

	HandleCudaError(cudaMemcpyToArray(gpOpacityArray, 0, 0, Opacity, TF_NO_SAMPLES * sizeof(float), cudaMemcpyHostToDevice));
	HandleCudaError(cudaBindTextureToArray(gTexOpacity, gpOpacityArray, ChannelDesc));
}

void UnbindTransferFunctionOpacity(void)
{
	HandleCudaError(cudaFreeArray(gpOpacityArray));
	gpOpacityArray = NULL;
	HandleCudaError(cudaUnbindTexture(gTexOpacity));
}

void BindTransferFunctionDiffuse(CTransferFunction& TransferFunctionDiffuse)
{
	gTexDiffuse.normalized		= true;
	gTexDiffuse.filterMode		= cudaFilterModeLinear;
	gTexDiffuse.addressMode[0]	= cudaAddressModeClamp;

	float4 Diffuse[TF_NO_SAMPLES];

	for (int i = 0; i < TF_NO_SAMPLES; i++)
	{
		Diffuse[i].x = TransferFunctionDiffuse.F((float)i * INV_TF_NO_SAMPLES).r;
		Diffuse[i].y = TransferFunctionDiffuse.F((float)i * INV_TF_NO_SAMPLES).g;
		Diffuse[i].z = TransferFunctionDiffuse.F((float)i * INV_TF_NO_SAMPLES).b;
	}

	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<float4>();
	
	if (gpDiffuseArray == NULL)
		HandleCudaError(cudaMallocArray(&gpDiffuseArray, &ChannelDesc, TF_NO_SAMPLES, 1));

	HandleCudaError(cudaMemcpyToArray(gpDiffuseArray, 0, 0, Diffuse, TF_NO_SAMPLES * sizeof(float4), cudaMemcpyHostToDevice));
	HandleCudaError(cudaBindTextureToArray(gTexDiffuse, gpDiffuseArray, ChannelDesc));
}

void UnbindTransferFunctionDiffuse(void)
{
	HandleCudaError(cudaFreeArray(gpDiffuseArray));
	gpDiffuseArray = NULL;
	HandleCudaError(cudaUnbindTexture(gTexDiffuse));
}

void BindTransferFunctionSpecular(CTransferFunction& TransferFunctionSpecular)
{
	gTexSpecular.normalized		= true;
	gTexSpecular.filterMode		= cudaFilterModeLinear;
	gTexSpecular.addressMode[0]	= cudaAddressModeClamp;

	float4 Specular[TF_NO_SAMPLES];

	for (int i = 0; i < TF_NO_SAMPLES; i++)
	{
		Specular[i].x = TransferFunctionSpecular.F((float)i * INV_TF_NO_SAMPLES).r;
		Specular[i].y = TransferFunctionSpecular.F((float)i * INV_TF_NO_SAMPLES).g;
		Specular[i].z = TransferFunctionSpecular.F((float)i * INV_TF_NO_SAMPLES).b;
	}

	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<float4>();
	
	if (gpSpecularArray == NULL)
		HandleCudaError(cudaMallocArray(&gpSpecularArray, &ChannelDesc, TF_NO_SAMPLES, 1));

	HandleCudaError(cudaMemcpyToArray(gpSpecularArray, 0, 0, Specular, TF_NO_SAMPLES * sizeof(float4), cudaMemcpyHostToDevice));
	HandleCudaError(cudaBindTextureToArray(gTexSpecular, gpSpecularArray, ChannelDesc));
}

void UnbindTransferFunctionSpecular(void)
{
	HandleCudaError(cudaFreeArray(gpSpecularArray));
	gpSpecularArray = NULL;
	HandleCudaError(cudaUnbindTexture(gTexSpecular));
}

void BindTransferFunctionRoughness(CTransferFunction& TransferFunctionRoughness)
{
	gTexRoughness.normalized		= true;
	gTexRoughness.filterMode		= cudaFilterModeLinear;
	gTexRoughness.addressMode[0]	= cudaAddressModeClamp;

	float Roughness[TF_NO_SAMPLES];

	for (int i = 0; i < TF_NO_SAMPLES; i++)
		Roughness[i] = TransferFunctionRoughness.F((float)i * INV_TF_NO_SAMPLES).r;
	
	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<float>();

	if (gpRoughnessArray == NULL)
		HandleCudaError(cudaMallocArray(&gpRoughnessArray, &ChannelDesc, TF_NO_SAMPLES, 1));

	HandleCudaError(cudaMemcpyToArray(gpRoughnessArray, 0, 0, Roughness, TF_NO_SAMPLES * sizeof(float),  cudaMemcpyHostToDevice));
	HandleCudaError(cudaBindTextureToArray(gTexRoughness, gpRoughnessArray, ChannelDesc));
}

void UnbindTransferFunctionRoughness(void)
{
	HandleCudaError(cudaFreeArray(gpRoughnessArray));
	gpRoughnessArray = NULL;
	HandleCudaError(cudaUnbindTexture(gTexRoughness));
}

void BindTransferFunctionEmission(CTransferFunction& TransferFunctionEmission)
{
	gTexEmission.normalized		= true;
	gTexEmission.filterMode		= cudaFilterModeLinear;
	gTexEmission.addressMode[0]	= cudaAddressModeClamp;

	float4 Emission[TF_NO_SAMPLES];

	for (int i = 0; i < TF_NO_SAMPLES; i++)
	{
		Emission[i].x = TransferFunctionEmission.F((float)i * INV_TF_NO_SAMPLES).r;
		Emission[i].y = TransferFunctionEmission.F((float)i * INV_TF_NO_SAMPLES).g;
		Emission[i].z = TransferFunctionEmission.F((float)i * INV_TF_NO_SAMPLES).b;
	}

	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<float4>();
	
	if (gpEmissionArray == NULL)
		HandleCudaError(cudaMallocArray(&gpEmissionArray, &ChannelDesc, TF_NO_SAMPLES, 1));

	HandleCudaError(cudaMemcpyToArray(gpEmissionArray, 0, 0, Emission, TF_NO_SAMPLES * sizeof(float4),  cudaMemcpyHostToDevice));
	HandleCudaError(cudaBindTextureToArray(gTexEmission, gpEmissionArray, ChannelDesc));
}

void UnbindTransferFunctionEmission(void)
{
	HandleCudaError(cudaFreeArray(gpEmissionArray));
	gpEmissionArray = NULL;
	HandleCudaError(cudaUnbindTexture(gTexEmission));
}

void BindConstants(CScene* pScene)
{
	const float3 AaBbMin = make_float3(pScene->m_BoundingBox.GetMinP().x, pScene->m_BoundingBox.GetMinP().y, pScene->m_BoundingBox.GetMinP().z);
	const float3 AaBbMax = make_float3(pScene->m_BoundingBox.GetMaxP().x, pScene->m_BoundingBox.GetMaxP().y, pScene->m_BoundingBox.GetMaxP().z);

	HandleCudaError(cudaMemcpyToSymbol(gAaBbMin, &AaBbMin, sizeof(float3)));
	HandleCudaError(cudaMemcpyToSymbol(gAaBbMax, &AaBbMax, sizeof(float3)));

	const float3 InvAaBbMin = make_float3(pScene->m_BoundingBox.GetInvMinP().x, pScene->m_BoundingBox.GetInvMinP().y, pScene->m_BoundingBox.GetInvMinP().z);
	const float3 InvAaBbMax = make_float3(pScene->m_BoundingBox.GetInvMaxP().x, pScene->m_BoundingBox.GetInvMaxP().y, pScene->m_BoundingBox.GetInvMaxP().z);

	HandleCudaError(cudaMemcpyToSymbol(gInvAaBbMin, &InvAaBbMin, sizeof(float3)));
	HandleCudaError(cudaMemcpyToSymbol(gInvAaBbMax, &InvAaBbMax, sizeof(float3)));

	const float IntensityMin		= pScene->m_IntensityRange.GetMin();
	const float IntensityMax		= pScene->m_IntensityRange.GetMax();
	const float IntensityRange		= pScene->m_IntensityRange.GetRange();
	const float IntensityInvRange	= 1.0f / IntensityRange;

	HandleCudaError(cudaMemcpyToSymbol(gIntensityMin, &IntensityMin, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gIntensityMax, &IntensityMax, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gIntensityRange, &IntensityRange, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gIntensityInvRange, &IntensityInvRange, sizeof(float)));

	const float StepSize		= pScene->m_StepSizeFactor * pScene->m_GradientDelta;
	const float StepSizeShadow	= pScene->m_StepSizeFactorShadow * pScene->m_GradientDelta;
	const float ScatteringHeadstart = pScene->m_ScatteringHeadstart * pScene->m_GradientDelta;

	//TODO remove
	//std::cout << "StepSize: " << StepSize << ", StepSizeShadow: " << StepSizeShadow << ", GradientDelta: "
	//	<< pScene->m_GradientDelta << ", DensityScale: " << pScene->m_DensityScale << ", GradientFactor: " << pScene->m_GradientFactor
	//	<< std::endl;
	/*CColorXyz tem = (CColorRgbHdr(1.0f).ToXYZ() * INV_4_PI_F) / INV_4_PI_F;
	CColorRgbHdr asdf;
	asdf.FromXYZ(tem.c[0], tem.c[1], tem.c[2]);*/
	//CColorRgbHdr asdf = CColorRgbHdr(INV_4_PI_F);
	//std::cout << asdf.r << ", " << asdf.g << ", " << asdf.b << std::endl;

	HandleCudaError(cudaMemcpyToSymbol(gStepSize, &StepSize, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gStepSizeShadow, &StepSizeShadow, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gScatteringHeadstart, &ScatteringHeadstart, sizeof(float)));


	const float DensityScale = pScene->m_DensityScale;

	HandleCudaError(cudaMemcpyToSymbol(gDensityScale, &DensityScale, sizeof(float)));
	
	const float GradientDelta		= 1.0f * pScene->m_GradientDelta;
	const float InvGradientDelta	= 1.0f / GradientDelta;
	const Vec3f GradientDeltaX(GradientDelta, 0.0f, 0.0f);
	const Vec3f GradientDeltaY(0.0f, GradientDelta, 0.0f);
	const Vec3f GradientDeltaZ(0.0f, 0.0f, GradientDelta);
	
	HandleCudaError(cudaMemcpyToSymbol(gGradientDelta, &GradientDelta, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gInvGradientDelta, &InvGradientDelta, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gGradientDeltaX, &GradientDeltaX, sizeof(Vec3f)));
	HandleCudaError(cudaMemcpyToSymbol(gGradientDeltaY, &GradientDeltaY, sizeof(Vec3f)));
	HandleCudaError(cudaMemcpyToSymbol(gGradientDeltaZ, &GradientDeltaZ, sizeof(Vec3f)));
	
	const int FilmWidth		= pScene->m_Camera.m_Film.GetWidth();
	const int Filmheight	= pScene->m_Camera.m_Film.GetHeight();
	const int FilmNoPixels	= pScene->m_Camera.m_Film.m_Resolution.GetNoElements();

	HandleCudaError(cudaMemcpyToSymbol(gFilmWidth, &FilmWidth, sizeof(int)));
	HandleCudaError(cudaMemcpyToSymbol(gFilmHeight, &Filmheight, sizeof(int)));
	HandleCudaError(cudaMemcpyToSymbol(gFilmNoPixels, &FilmNoPixels, sizeof(int)));

	const int FilterWidth = 1;

	HandleCudaError(cudaMemcpyToSymbol(gFilterWidth, &FilterWidth, sizeof(int)));

	const float FilterWeights[10] = { 0.11411459588254977f, 0.08176668094332218f, 0.03008028089187349f, 0.01f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };

	HandleCudaError(cudaMemcpyToSymbol(gFilterWeights, &FilterWeights, 10 * sizeof(float)));

	const float Gamma		= pScene->m_Camera.m_Film.m_Gamma;
	const float InvGamma	= 1.0f / Gamma;
	const float Exposure	= pScene->m_Camera.m_Film.m_Exposure;
	const float InvExposure	= 1.0f / Exposure;
	// TODO remove: std::cout << "invExposure: " << InvExposure << std::endl;
	HandleCudaError(cudaMemcpyToSymbol(gExposure, &Exposure, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gInvExposure, &InvExposure, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gGamma, &Gamma, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gInvGamma, &InvGamma, sizeof(float)));

	HandleCudaError(cudaMemcpyToSymbol(gDenoiseEnabled, &pScene->m_DenoiseParams.m_Enabled, sizeof(bool)));
	HandleCudaError(cudaMemcpyToSymbol(gDenoiseWindowRadius, &pScene->m_DenoiseParams.m_WindowRadius, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gDenoiseInvWindowArea, &pScene->m_DenoiseParams.m_InvWindowArea, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gDenoiseNoise, &pScene->m_DenoiseParams.m_Noise, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gDenoiseWeightThreshold, &pScene->m_DenoiseParams.m_WeightThreshold, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gDenoiseLerpThreshold, &pScene->m_DenoiseParams.m_LerpThreshold, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gDenoiseLerpC, &pScene->m_DenoiseParams.m_LerpC, sizeof(float)));

	const float NoIterations	= pScene->GetNoIterations();
	const float InvNoIterations = 1.0f / __max(1.0f, NoIterations);

	HandleCudaError(cudaMemcpyToSymbol(gNoIterations, &NoIterations, sizeof(float)));
	HandleCudaError(cudaMemcpyToSymbol(gInvNoIterations, &InvNoIterations, sizeof(float)));
}

void Render(CScene& Scene, CTiming& RenderImage, CTiming& BlurImage, CTiming& PostProcessImage, CTiming& DenoiseImage, curandState* pStates)
{
	CScene* pDevScene = NULL;

	HandleCudaError(cudaMalloc(&pDevScene, sizeof(CScene)));
	HandleCudaError(cudaMemcpy(pDevScene, &Scene, sizeof(CScene), cudaMemcpyHostToDevice));

	if (Scene.m_Camera.m_Focus.m_Type == 0)
		Scene.m_Camera.m_Focus.m_FocalDistance = NearestIntersection(pDevScene);

	HandleCudaError(cudaMemcpy(pDevScene, &Scene, sizeof(CScene), cudaMemcpyHostToDevice));

	CCudaView* pDevView = NULL;

	HandleCudaError(cudaMalloc(&pDevView, sizeof(CCudaView)));
	HandleCudaError(cudaMemcpy(pDevView, &gRenderCanvasView, sizeof(CCudaView), cudaMemcpyHostToDevice));

	//std::cout << INV_4_PI_F << std::endl;
	
	CCudaTimer TmrRender;
	
	switch (Scene.m_AlgorithmType)
	{
		case -1: {
			GiveGreen(&Scene, pDevScene, pDevView);
			break;
		}
		case 0:
		{
			SingleScattering(&Scene, pDevScene, pDevView);
			break;
		}

		case 5: {
			MultipleScattering(&Scene, pDevScene, pDevView);
			break;
		}

		case 6: {
			MultipleScattering(&Scene, pDevScene, pDevView);
			break;
		}

		case 7: {
			MultipleScattering(&Scene, pDevScene, pDevView);
			break;
		}

		case 1:
		{
			//MultipleScattering(&Scene, pDevScene);
			MultipleScattering(&Scene, pDevScene, pDevView);
			break;
		}

		case 2:
		{
			//int SizeResults = Scene.m_Camera.m_Film.m_Resolution.GetResX() * Scene.m_Camera.m_Film.m_Resolution.GetResY() * sizeof(Vec3f);
			//int SizeResults = m_nrPoints * sizeof(Vec3f);
			Vec3f* pDevResults = NULL;
			//HandleCudaError(cudaMalloc(&pDevResults, SizeResults));


			PreCalculatedScattering(&Scene, pDevScene, pDevView, m_DevPoints, m_nrPoints, m_DevConnections, m_DevColours , pDevResults);
			
			/*size_t free = -1;
			size_t total = -1;
			HandleCudaError(cudaMemGetInfo(&free, &total));
			cout << (free / (sizeof(Vec3f) + sizeof(CColorXyza))) << ", " << total << endl;*/

			/*Vec3f* pResults;
			pResults = (Vec3f *)malloc(SizeResults);
			HandleCudaError(cudaMemcpy(pResults, pDevResults, SizeResults, cudaMemcpyDeviceToHost));
			
			int xWidth = Scene.m_Camera.m_Film.m_Resolution.GetResX();
			for (int i = 0; i < Scene.m_Camera.m_Film.m_Resolution.GetResX() * Scene.m_Camera.m_Film.m_Resolution.GetResY(); i++) {
				int x = xWidth != 0 ? i % xWidth : -1;
				int y = xWidth != 0 ? i / xWidth : -1;
				std::cout << "i: " << i << ", X: " << x << ", Y: " << y << ", p: " << pResults[i].x << ", " << pResults[i].y << ", " << pResults[i].z << endl;
			}
			char x;
			std::cin >> x;

			free(pResults);*/
			//HandleCudaError(cudaFree(pDevResults));
			




/*int SizePoints = Scene.m_Camera.m_Film.m_Resolution.GetResX() * Scene.m_Camera.m_Film.m_Resolution.GetResY() * sizeof(Vec3f);
Vec3f* pDevPoints;
HandleCudaError(cudaMalloc(&pDevPoints, SizePoints));

getScatterInfo(&Scene, pDevScene, pDevView, pDevPoints, pStates);

Vec3f* pPoints;
pPoints = (Vec3f *)malloc(SizePoints);
HandleCudaError(cudaMemcpy(pPoints, pDevPoints, SizePoints, cudaMemcpyDeviceToHost));
*/
//int nrBins = 100;
//float binSize = 1.f / nrBins;
//int* bins;
//bins = (int*)malloc(nrBins * sizeof(int));
//for (int i = 0; i < nrBins; i++)
//	bins[i] = 0;
//int samples = 0;
//float min = 10.f;
//float max = -10.f;

////for (int i = 0; i < Scene.m_Camera.m_Film.m_Resolution.GetResX() * Scene.m_Camera.m_Film.m_Resolution.GetResY(); i++) {
//for (int i = 0; i <m_nrPoints; i++) {
//	//std::cout << "i: " << i << ", p: " << pResults[i].x << ", " << pResults[i].y << ", " << pResults[i].z << endl;

//	float sample = pResults[i].x;
//	int binId = sample / binSize;
//	if (binId < nrBins) {
//		if (binId < 0 || binId >= nrBins)
//			continue;
//		bins[binId] = bins[binId] + 1;
//		samples++;
//		min = Fminf(sample, min);
//		max = Fmaxf(sample, max);
//	}
//}

//for (int i = 0; i < nrBins; i++) {
//	std::cout << (i * binSize) << ": " << bins[i] << endl;
//}

//std::cout << "Samples: " << samples << endl;
//std::cout << "Min: " << min << endl;
//std::cout << "Max: " << max << endl;
//
//

//char x;
//std::cin >> x;

break;
		}

		case 3: {
			SingleScattering(&Scene, pDevScene, pDevView);
			break;
		}
	}

	RenderImage.AddDuration(TmrRender.ElapsedTime());

	CCudaTimer TmrBlur;
	if (Scene.m_PostProcessingSteps & PostProcessingStepsEnum::BLUR) {
		Blur(&Scene, pDevScene, pDevView);
	}
	BlurImage.AddDuration(TmrBlur.ElapsedTime());

	CCudaTimer TmrPostProcess;
	if (Scene.m_PostProcessingSteps & PostProcessingStepsEnum::ESTIMATE) {
		Estimate(&Scene, pDevScene, pDevView);
	}
	else {
		EstimateCopy(&Scene, pDevScene, pDevView);
	}
	PostProcessImage.AddDuration(TmrPostProcess.ElapsedTime());

	if (Scene.m_PostProcessingSteps & PostProcessingStepsEnum::TONE_MAP) {
		ToneMap(&Scene, pDevScene, pDevView);
	}
	else {
		ToneMapCopy(&Scene, pDevScene, pDevView);
	}

	CCudaTimer TmrDenoise;
	if (Scene.m_PostProcessingSteps & PostProcessingStepsEnum::DENOISE) {
		Denoise(&Scene, pDevScene, pDevView);
	}
	else {
		//pDevView->m_DisplayEstimateRgbLdr = pDevView->m_EstimateRgbaLdr;
		DenoiseCopy(&Scene, pDevScene, pDevView);
	}
	DenoiseImage.AddDuration(TmrDenoise.ElapsedTime());

	HandleCudaError(cudaFree(pDevScene));
	HandleCudaError(cudaFree(pDevView));
}

void InitOpacityGradient(CScene& Scene) {
	CCudaTimer TmrinitGradient;

	cudaExtent Extent;
	Extent.width = Scene.m_Resolution[0];
	Extent.height = Scene.m_Resolution[1];
	Extent.depth = Scene.m_Resolution[2];

	int sizeOpacityGradient = Scene.m_Resolution.GetNoElements() * sizeof(float4);
	int sizeOpacityGradientMagnitude = Scene.m_Resolution.GetNoElements() * sizeof(float);

	float4* pDevOpacityGradient;
	HandleCudaError(cudaMalloc(&pDevOpacityGradient, sizeOpacityGradient));

	float* pDevOpacityGradientMagnitude1D;
	HandleCudaError(cudaMalloc(&pDevOpacityGradientMagnitude1D, sizeOpacityGradientMagnitude));

	GetOpacityGradientTexture(Scene.m_Resolution, pDevOpacityGradient, pDevOpacityGradientMagnitude1D);
	
	float4* pOpacityGradient1D;
	pOpacityGradient1D = (float4 *)malloc(sizeOpacityGradient);
	HandleCudaError(cudaMemcpy(pOpacityGradient1D, pDevOpacityGradient, sizeOpacityGradient, cudaMemcpyDeviceToHost));

	float* pOpacityGradientMagnitude1D;
	pOpacityGradientMagnitude1D = (float*)malloc(sizeOpacityGradientMagnitude);
	HandleCudaError(cudaMemcpy(pOpacityGradientMagnitude1D, pDevOpacityGradientMagnitude1D, sizeOpacityGradientMagnitude, cudaMemcpyDeviceToHost));

	// Normalize the magnitudes
	/*float max = -1;
	float min = std::numeric_limits<float>::max();
	for (int i = 0; i < Scene.m_Resolution.GetNoElements(); i++) {
		max = Fmaxf(max, pOpacityGradientMagnitude1D[i]);
		min = Fminf(min, pOpacityGradientMagnitude1D[i]);
	}
	for (int i = 0; i < Scene.m_Resolution.GetNoElements(); i++) {
		pOpacityGradientMagnitude1D[i] = pOpacityGradientMagnitude1D[i] / max;
	}*/
	
	// Transfer the 1D Gradient array to a 3D texture
	cudaChannelFormatDesc ChannelDesc = cudaCreateChannelDesc<float4>();
	cudaArray* pCuOpacityGradient = NULL;
	HandleCudaError(cudaMalloc3DArray(&pCuOpacityGradient, &ChannelDesc, Extent));

	cudaMemcpy3DParms CopyParamsGradient = { 0 };

	CopyParamsGradient.srcPtr = make_cudaPitchedPtr(pOpacityGradient1D, Extent.width * sizeof(float4), Extent.width, Extent.height);
	CopyParamsGradient.dstArray = pCuOpacityGradient;
	CopyParamsGradient.extent = Extent;
	CopyParamsGradient.kind = cudaMemcpyHostToDevice;

	HandleCudaError(cudaMemcpy3D(&CopyParamsGradient));

	gTexOpacityGradient.normalized = false;
	gTexOpacityGradient.filterMode = cudaFilterModeLinear;
	gTexOpacityGradient.addressMode[0] = cudaAddressModeClamp;
	gTexOpacityGradient.addressMode[1] = cudaAddressModeClamp;
	gTexOpacityGradient.addressMode[2] = cudaAddressModeClamp;

	HandleCudaError(cudaBindTextureToArray(gTexOpacityGradient, pCuOpacityGradient, ChannelDesc));


	// Bind the 1D magnitude array to 3D texture
	cudaChannelFormatDesc ChannelDescMag = cudaCreateChannelDesc<float>();
	cudaArray* pCuOpacityMagnitudeNormalized = NULL;
	HandleCudaError(cudaMalloc3DArray(&pCuOpacityMagnitudeNormalized, &ChannelDescMag, Extent));

	cudaMemcpy3DParms CopyParams = { 0 };

	CopyParams.srcPtr = make_cudaPitchedPtr(pOpacityGradient1D, Extent.width * sizeof(float), Extent.width, Extent.height);
	CopyParams.dstArray = pCuOpacityMagnitudeNormalized;
	CopyParams.extent = Extent;
	CopyParams.kind = cudaMemcpyHostToDevice;

	HandleCudaError(cudaMemcpy3D(&CopyParams));

	gTexOpacityMagnitudeNormalized.normalized = true;
	gTexOpacityMagnitudeNormalized.filterMode = cudaFilterModeLinear;
	gTexOpacityMagnitudeNormalized.addressMode[0] = cudaAddressModeClamp;
	gTexOpacityMagnitudeNormalized.addressMode[1] = cudaAddressModeClamp;
	gTexOpacityMagnitudeNormalized.addressMode[2] = cudaAddressModeClamp;

	HandleCudaError(cudaBindTextureToArray(gTexOpacityMagnitudeNormalized, pCuOpacityMagnitudeNormalized, ChannelDescMag));

	/*gTexDensity = gTexOpacityMagnitudeNormalized;

	// Make 1to1 opacity texture
	texture<float, cudaTextureType1D, cudaReadModeElementType> gTexOpacityDummy;
	gTexOpacityDummy.normalized = true;
	gTexOpacityDummy.filterMode = cudaFilterModeLinear;
	gTexOpacityDummy.addressMode[0] = cudaAddressModeClamp;

	float Opacity[2];
	Opacity[0] = 0;
	Opacity[1] = 1;

	cudaChannelFormatDesc ChannelDescDummy = cudaCreateChannelDesc<float>();

	cudaArray* gpOpacityArrayDummy = NULL;
	if (gpOpacityArrayDummy == NULL)
		HandleCudaError(cudaMallocArray(&gpOpacityArrayDummy, &ChannelDescDummy, 2, 1));

	HandleCudaError(cudaMemcpyToArray(gpOpacityArrayDummy, 0, 0, Opacity, TF_NO_SAMPLES * sizeof(float), cudaMemcpyHostToDevice));
	HandleCudaError(cudaBindTextureToArray(gTexOpacityDummy, gpOpacityArrayDummy, ChannelDescDummy));*/

	cout << "Getting opacity gradient magnitude: " << TmrinitGradient.ElapsedTime() << " ms" << endl;

	int nrPoints = 150000;

	Vec3i* points = getPointsOpacityGradientMagnitudeBased(pOpacityGradientMagnitude1D, Scene.m_Resolution, nrPoints);

	int maxX = -1;
	int maxY = -1;
	int maxZ = -1;
	for (int i = 0; i < nrPoints; i++) {
		maxX = max(maxX, points[i].x);
		maxY = max(maxY, points[i].y);
		maxZ = max(maxZ, points[i].z);
	}
	cout << "Max(x,y,z)" << maxX << ", " << maxY << ", " << maxZ << endl;
	cout << "Res(x,y,z)" << Scene.m_Resolution.GetResX() << ", " << Scene.m_Resolution.GetResY() << ", " << Scene.m_Resolution.GetResZ() << endl;

	// override density buffer;
	overridDensity(points, Scene.m_Resolution, nrPoints);

	// create graph

}

void overridDensity(Vec3i* points, CResolution3D resolution, int nrPoits) {
	CCudaTimer TmrReplaceTexture;

	short* pFakeDensityBuffer;
	pFakeDensityBuffer = (short*)malloc(resolution.GetNoElements() * sizeof(short));

	for (int i = 0; i < resolution.GetNoElements(); i++) {
		pFakeDensityBuffer[i] = -1000;
	}

	for (int i = 0; i < nrPoits; i++) {
		int index = Index3To1(points[i], resolution);
		pFakeDensityBuffer[index] = 3095;
	}

	UnbindDensityBuffer();

	cudaExtent Extent;
	Extent.width = resolution[0];
	Extent.height = resolution[1];
	Extent.depth = resolution[2];

	BindDensityBuffer(pFakeDensityBuffer, Extent);

	cout << "Replacing Texture: " << TmrReplaceTexture.ElapsedTime() << " ms" << endl;
}

Vec3i* getPointsOpacityGradientMagnitudeBased(float* opacityGradientMagnitudes, CResolution3D resolution, int nrPoints) {
	CCudaTimer TmrMakePdf;

	// Create PDF's
	vector<float> pdfX(resolution.GetResX());
	vector<vector<float>> pdfY(resolution.GetResX(), vector<float>(resolution.GetResY()));
	vector<vector<vector<float>>> pdfZ(resolution.GetResX(), vector<vector<float>>(resolution.GetResY(), vector<float>(resolution.GetResZ())));

	for (int x = 0; x < resolution.GetResX(); x++) {
		vector<float> pdfYcur = pdfY.at(x);
		for (int y = 0; y < resolution.GetResY(); y++) {
			vector<float> pdfZcur = pdfZ.at(x).at(y);
			for (int z = 0; z < resolution.GetResZ(); z++) {
				int index = Index3To1(x, y, z, resolution);
				float val = opacityGradientMagnitudes[index];
				
				if (y == 0 && z == 0 && x != 0)
					pdfX.at(x) = pdfX.at(x - 1) + val;
				else
					pdfX.at(x) += val;
				
				if (z == 0 && y != 0)
					pdfYcur.at(y) = pdfYcur.at(y - 1) + val;
				else
					pdfYcur.at(y) += val;
				
				if (z != 0)
					pdfZcur.at(z) = pdfZcur.at(z - 1) + val;
				else
					pdfZcur.at(z) = val;
			}
			pdfZ.at(x).at(y) = pdfZcur;
		}
		pdfY.at(x) = pdfYcur;
	}

	cout << "Making PDF's: " << TmrMakePdf.ElapsedTime() << " ms" << endl;

	/*cout << fixed;
	for (int i = 0; i < resolution.GetResX(); i++) {
		cout << pdfX[i] << endl;
	}*/

	//int nrPoints = 150;

	CCudaTimer TmrGetPoints;

	Vec3i* points;
	points = (Vec3i*)malloc(nrPoints * sizeof(Vec3i));

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0, 1.0);
	for (int n = 0; n < nrPoints; n++) {
		double val = dis(gen) * pdfX.back();

		int x, y, z;

		for (x = 0; x < resolution.GetResX(); x++) {
			if (pdfX[x] > val) {
				break;
			}
		}

		vector<float> pdfYcur = pdfY.at(x);
		val = dis(gen) * pdfYcur.back();

		for (y = 0; y < resolution.GetResY(); y++) {
			if (pdfYcur[y] > val) {
				break;
			}
		}

		vector<float> pdfZcur = pdfZ.at(x).at(y);
		val = dis(gen) * pdfZcur.back();

		for (z = 0; z < resolution.GetResZ(); z++) {
			if (pdfZcur[z] > val) {
				break;
			}
		}

		points[n] = Vec3i(x, y, z);
	}
	

	/*for (int i = 0; i < nrPoints; i++) {
		Vec3i p = points[i];
		cout << p.x << ", " << p.y << ", " << p.z << endl;
	}*/
	cout << "Generating points: " << TmrGetPoints.ElapsedTime() << endl;

	return points;
}


void InitPreCalculatedCore(CScene& Scene, short* pDensityBuffer) {
	if (m_DevPoints != NULL)
		HandleCudaError(cudaFree(m_DevPoints));
	if (m_DevColours != NULL)
		HandleCudaError(cudaFree(m_DevColours));
	if (m_DevConnections != NULL)
		HandleCudaError(cudaFree(m_DevConnections));

	CScene* pDevScene = NULL;
	HandleCudaError(cudaMalloc(&pDevScene, sizeof(CScene)));
	HandleCudaError(cudaMemcpy(pDevScene, &Scene, sizeof(CScene), cudaMemcpyHostToDevice));

	int nrPoints = Scene.m_Resolution.GetNoElements() / 500;
	//int nrPoints = Scene.m_Camera.m_Film.GetWidth() * Scene.m_Camera.m_Film.GetHeight();
	int PointPerState = 10;
	int nrThreads = (nrPoints + (PointPerState - 1)) / PointPerState;

	int SizePoints = nrPoints * sizeof(Vec3f);

	cout << "Total Points: " << Scene.m_Resolution.GetNoElements() << ", Generated Points: " << nrPoints << endl;

	Vec3f* pDevPoints;
	HandleCudaError(cudaMalloc(&pDevPoints, SizePoints));

	int StateSize = nrThreads * sizeof(curandState);

	curandState* pDevStates;
	HandleCudaError(cudaMalloc(&pDevStates, StateSize));

	int SeedSize = nrThreads * sizeof(int);

	int* Seeds;
	Seeds = (int *)malloc(SeedSize);
	random_ints(Seeds, nrThreads);

	int* pDevSeeds;
	HandleCudaError(cudaMalloc(&pDevSeeds, SeedSize));
	HandleCudaError(cudaMemcpy(pDevSeeds, Seeds, SeedSize, cudaMemcpyHostToDevice));

	SetupStates(pDevStates, pDevSeeds, nrThreads);

	GeneratePoints(pDevScene, pDevStates, pDevPoints, nrThreads, PointPerState, nrPoints);

	Vec3f* pPoints;
	pPoints = (Vec3f *)malloc(SizePoints);
	HandleCudaError(cudaMemcpy(pPoints, pDevPoints, SizePoints, cudaMemcpyDeviceToHost));

	
	HandleCudaError(cudaFree(pDevStates));
	HandleCudaError(cudaFree(pDevSeeds));
	free(Seeds);

	// generate colours/connections

	int seedSizeRNG = nrPoints * sizeof(unsigned int);
	unsigned int* pSeedsCRNG = (unsigned int*)malloc(seedSizeRNG);

	// Set  seeds1
	for (int i = 0; i < nrPoints; i++)
		pSeedsCRNG[i] = rand();

	unsigned int* pDevSeeds1CRNG;
	HandleCudaError(cudaMalloc(&pDevSeeds1CRNG, seedSizeRNG));
	HandleCudaError(cudaMemcpy(pDevSeeds1CRNG, pSeedsCRNG, seedSizeRNG, cudaMemcpyHostToDevice));

	// Set seeds2
	for (int i = 0; i < nrPoints; i++)
		pSeedsCRNG[i] = rand();

	unsigned int* pDevSeeds2CRNG;
	HandleCudaError(cudaMalloc(&pDevSeeds2CRNG, seedSizeRNG));
	HandleCudaError(cudaMemcpy(pDevSeeds2CRNG, pSeedsCRNG, seedSizeRNG, cudaMemcpyHostToDevice));

	free(pSeedsCRNG);

	float* pDevConnections = NULL;
	//HandleCudaError(cudaMalloc(&pDevConnections, nrPoints * nrPoints * sizeof(float)));

	CColorXyza* pDevPointColours;
	HandleCudaError(cudaMalloc(&pDevPointColours, nrPoints * sizeof(CColorXyza)));


	makeConnections(pDevScene, pDevPoints, nrPoints, pDevSeeds1CRNG, pDevSeeds2CRNG, pDevConnections, pDevPointColours);

	/*CColorXyza* pPointColours;
	pPointColours = (CColorXyza*)malloc(nrPoints * sizeof(CColorXyza));
	HandleCudaError(cudaMemcpy(pPointColours, pDevPointColours, nrPoints * sizeof(CColorXyza), cudaMemcpyDeviceToHost));*/

	/*float* pConnections;
	pConnections = (float*)malloc(nrPoints * nrPoints * sizeof(float));
	HandleCudaError(cudaMemcpy(pConnections, pDevConnections, nrPoints * nrPoints * sizeof(float), cudaMemcpyDeviceToHost));*/


	// Store data
	m_nrPoints = nrPoints;
	m_DevPoints = pDevPoints;
	m_DevColours = pDevPointColours;
	m_DevConnections = pDevConnections;

	HandleCudaError(cudaFree(pDevSeeds1CRNG));
	HandleCudaError(cudaFree(pDevSeeds2CRNG));
	
	HandleCudaError(cudaFree(pDevScene));
}

/// <summary>
/// Initializes N curandstates
/// </summary>
/// <param name="N">Number of required states</param>
curandState* InitStates(int N) {
	int StateSize = N * sizeof(curandState);

	curandState* pDevStates;
	HandleCudaError(cudaMalloc(&pDevStates, StateSize));

	int SeedSize = N * sizeof(int);
	
	int* pSeeds;
	pSeeds = (int *)malloc(SeedSize);
	random_ints(pSeeds, N);

	int* pDevSeeds;
	HandleCudaError(cudaMalloc(&pDevSeeds, SeedSize));
	HandleCudaError(cudaMemcpy(pDevSeeds, pSeeds, SeedSize, cudaMemcpyHostToDevice));

	//SetupStates(pDevStates, pDevSeeds, N);
	SetupStatesSingleSeed(pDevStates, rand() , N);

	HandleCudaError(cudaFree(pDevSeeds));
	free(pSeeds);

	return pDevStates;
}

void random_ints(int* a, int N)
{
	int i;
	for (i = 0; i < N; ++i)
		a[i] = rand();
}

void CreateIlluminanceTextureCore(CScene& Scene, float* pIlluminanceBuffer) {
	
	CResolution3D resolution = Scene.m_Resolution;

	/*CCudaRandomBuffer2D* pDevRandom1 = NULL;
	CCudaRandomBuffer2D random1;
	random1.Resize(CResolution2D(resolution.GetResX(), resolution.GetResY()));

	HandleCudaError(cudaMalloc(&pDevRandom1, sizeof(CCudaRandomBuffer2D)));
	HandleCudaError(cudaMemcpy(pDevRandom1, &random1, sizeof(CCudaRandomBuffer2D), cudaMemcpyHostToDevice));
	random1.Free();

	CCudaRandomBuffer2D* pDevRandom2 = NULL;
	CCudaRandomBuffer2D random2;
	random2.Resize(CResolution2D(resolution.GetResX(), resolution.GetResY()));

	HandleCudaError(cudaMalloc(&pDevRandom2, sizeof(CCudaRandomBuffer2D)));
	HandleCudaError(cudaMemcpy(pDevRandom2, &random2, sizeof(CCudaRandomBuffer2D), cudaMemcpyHostToDevice));
	random2.Free();*/

	int sizeSeeds = resolution.GetNoElements() * sizeof(unsigned int);
	unsigned int* seeds1;
	seeds1 = (unsigned int *)malloc(sizeSeeds);
	for (int i = 0; i < resolution.GetNoElements(); i++) {
		seeds1[i] = rand();
	}
	unsigned int* pDevSeeds1;
	HandleCudaError(cudaMalloc(&pDevSeeds1, sizeSeeds));
	HandleCudaError(cudaMemcpy(pDevSeeds1, seeds1, sizeSeeds, cudaMemcpyHostToDevice));

	unsigned int* seeds2;
	seeds2 = (unsigned int *)malloc(sizeSeeds);
	for (int i = 0; i < resolution.GetNoElements(); i++) {
		seeds2[i] = rand();
	}
	unsigned int* pDevSeeds2;
	HandleCudaError(cudaMalloc(&pDevSeeds2, sizeSeeds));
	HandleCudaError(cudaMemcpy(pDevSeeds2, seeds1, sizeSeeds, cudaMemcpyHostToDevice));

	CScene* pDevScene = NULL;
	HandleCudaError(cudaMalloc(&pDevScene, sizeof(CScene)));
	HandleCudaError(cudaMemcpy(pDevScene, &Scene, sizeof(CScene), cudaMemcpyHostToDevice));

	int SizeIlluminance = resolution.GetNoElements() * sizeof(float);
	float* pDevIlluminanceTexture;
	HandleCudaError(cudaMalloc(&pDevIlluminanceTexture, SizeIlluminance));

	for (int i = 0; i < 1; i++) {
		CCudaTimer timerIlluminanceRound;
		std::cout << "Working on iteration: " << i << std::endl;
		//GetIlluminationTexture(pDevScene, resolution, pDevIlluminanceTexture, pDevRandom1, pDevRandom2);
		GetIlluminationTexture(pDevScene, resolution, pDevIlluminanceTexture, pDevSeeds1, pDevSeeds2);


		std::cout << "Pass time: " << timerIlluminanceRound.ElapsedTime() << " ms" << std::endl;
	}

	HandleCudaError(cudaMemcpy(pIlluminanceBuffer, pDevIlluminanceTexture, SizeIlluminance, cudaMemcpyDeviceToHost));

	//free memory
	//HandleCudaError(cudaFree(pDevRandom1));
	//HandleCudaError(cudaFree(pDevRandom2));
	HandleCudaError(cudaFree(pDevSeeds1));
	free(seeds1);
	HandleCudaError(cudaFree(pDevSeeds2));
	free(seeds2);
	HandleCudaError(cudaFree(pDevScene));
	HandleCudaError(cudaFree(pDevIlluminanceTexture));
}