/*
	Copyright (c) 2011, T. Kroes <t.kroes@tudelft.nl>
	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

	- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
	- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
	- Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
	
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "Stable.h"

#include "RenderThread.h"
#include "CudaUtilities.h"
#include "MainWindow.h"
#include "LoadSettingsDialog.h"
#include "Lighting.h"
#include "Timing.h"

// std
#include <string>
#include <iostream>
#include <queue>
#include <chrono>

// CUDA kernels
#include "Core.cuh"

// VTK
#include <vtkSmartPointer.h>
#include <vtkMetaImageReader.h>
#include <vtkImageCast.h>
#include <vtkImageResample.h>
#include <vtkImageData.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkCallbackCommand.h>
#include <vtkImageAccumulate.h>
#include <vtkIntArray.h>
#include <vtkImageShiftScale.h>
#include <vtkErrorCode.h>
#include <vtkImageGradient.h>
#include <vtkExtractVectorComponents.h>
#include <vtkMetaImageWriter.h>

// Render thread
QRenderThread* gpRenderThread = NULL;
QFrameBuffer gFrameBuffer;

QMutex gSceneMutex;

int gCurrentDeviceID = 0;

curandState* pDevStates = nullptr;

QFrameBuffer::QFrameBuffer(void) :
	m_pPixels(NULL),
	m_Width(0),
	m_Height(0),
	m_NoPixels(0),
	m_Mutex()
{
}

QFrameBuffer::QFrameBuffer(const QFrameBuffer& Other)
{
	*this = Other;
}

QFrameBuffer& QFrameBuffer::operator=(const QFrameBuffer& Other)
{
	const bool Dirty = m_Width != Other.m_Width || m_Height != Other.m_Height;

	m_Width		= Other.m_Width;
	m_Height	= Other.m_Height;
	m_NoPixels	= Other.m_NoPixels;

	if (Other.m_pPixels != NULL)
	{
		const int Size = 3 * m_NoPixels * sizeof(unsigned char);

		if (Dirty)
		{
			free(m_pPixels);
			m_pPixels = (unsigned char*)malloc(Size);
		}

		memcpy(m_pPixels, Other.m_pPixels, Size); 
	}
	else
	{
		m_pPixels = NULL;
	}

	return *this;
}

QFrameBuffer::~QFrameBuffer(void)
{
	free(m_pPixels);
}

void QFrameBuffer::Set(unsigned char* pPixels, const int& Width, const int& Height)
{
	const bool Dirty = Width != m_Width || Height != m_Height;

	m_Width		= Width;
	m_Height	= Height;
	m_NoPixels	= m_Width * m_Height;

	if (m_NoPixels <= 0)
		return;

	const int Size = 3 * m_NoPixels * sizeof(unsigned char);

	if (Dirty)
	{
		free(m_pPixels);
		m_pPixels = (unsigned char*)malloc(Size);
	}

	memcpy(m_pPixels, pPixels, Size);
}

QRenderThread::QRenderThread(const QString& FileName, QObject* pParent /*= NULL*/) :
	QThread(pParent),
	m_FileName(FileName),
	m_pRenderImage(NULL),
	m_pDensityBuffer(NULL),
	m_pGradientMagnitudeBuffer(NULL),
	m_Abort(false),
	m_Pause(false),
	m_SaveFrames(),
	m_SaveBaseName("phase_function"),
	m_doTests(false),
	m_startTesting(false)
{
//	m_SaveFrames << 0 << 100 << 200;
}

QRenderThread::QRenderThread(const QRenderThread& Other)
{
	*this = Other;
}

QRenderThread::~QRenderThread(void)
{
	free(m_pDensityBuffer);
}

QRenderThread& QRenderThread::operator=(const QRenderThread& Other)
{
	m_FileName					= Other.m_FileName;
	m_pRenderImage				= Other.m_pRenderImage;
	m_pDensityBuffer			= Other.m_pDensityBuffer;
	m_pGradientMagnitudeBuffer	= Other.m_pGradientMagnitudeBuffer;
	m_Abort						= Other.m_Abort;
	m_Pause						= Other.m_Pause;
	m_SaveFrames				= Other.m_SaveFrames;
	m_SaveBaseName				= Other.m_SaveBaseName;

	return *this;
}

void QRenderThread::StartTesting(QString directory) {
	if (!isTesting()) {
		PauseRendering(true);

		m_startTesting = true;
		m_TestDir = directory;
		m_SaveFrames.append({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 });
		QDir dir = QDir(m_TestDir);
		dir.mkdir("images");

		QFile measurements(m_TestDir + "/measurements.csv");
		if (measurements.open(QIODevice::WriteOnly | QIODevice::Text)) {
			QTextStream stream(&measurements);
			QStringList data = {
				"iterations",
				"render_time",
				"blur_time",
				"postprocessing_time",
				"denoise_time",
				"used_memory_mb"
			};
			stream << data.join(",");
			measurements.close();
		}

		Log("Saving tests to: " + m_TestDir, "control");

		PauseRendering(false);
	}
	else {
		Log("Already performing a test", "control");
	}
}

void QRenderThread::run()
{
	if (!SetCudaDevice(gCurrentDeviceID))
		return;

	ResetDevice();

	CScene SceneCopy;
	
 	gScene.m_Camera.m_SceneBoundingBox = gScene.m_BoundingBox;
	gScene.m_Camera.SetViewMode(ViewModeFront);
 	gScene.m_Camera.Update();

	// Force the render thread to allocate the necessary buffers, do not remove this line
	gScene.m_DirtyFlags.SetFlag(FilmResolutionDirty | CameraDirty);

	gStatus.SetStatisticChanged("Memory", "CUDA Memory", "", "", "memory");
	gStatus.SetStatisticChanged("Memory", "Host Memory", "", "", "memory");

	cudaExtent Res;
	Res.width = gScene.m_Resolution[0];
	Res.height = gScene.m_Resolution[1];
	Res.depth = gScene.m_Resolution[2];

	// Bind density buffer to texture
	Log("Copying density volume to device", "grid");
	gStatus.SetStatisticChanged("CUDA Memory", "Density Buffer", QString::number(gScene.m_Resolution.GetNoElements() * sizeof(short) / MB, 'f', 2), "MB");
	BindDensityBuffer((short*)m_pDensityBuffer, Res, cudaFilterModeLinear);
	m_CurFilterMode = cudaFilterModeLinear;

	// Bind gradient magnitude buffer to texture
	Log("Copying gradient magnitude to device", "grid");
	gStatus.SetStatisticChanged("CUDA Memory", "Gradient Magnitude Buffer", QString::number(gScene.m_Resolution.GetNoElements() * sizeof(short) / MB, 'f', 2), "MB");
	BindGradientMagnitudeBuffer((short*)m_pGradientMagnitudeBuffer, Res);

	gStatus.SetStatisticChanged("Performance", "Timings", "");
	gStatus.SetStatisticChanged("Camera", "", "");
	gStatus.SetStatisticChanged("CUDA Memory", "Scene", QString::number(sizeof(CScene) / MB, 'f', 2), "MB");
	gStatus.SetStatisticChanged("CUDA Memory", "Frame Buffers", "", "");

	// Let others know that we are starting with rendering
	gStatus.SetRenderBegin();
	
	Log("Device memory: " + QString::number(GetUsedCudaMemory() / MB, 'f', 2) + "/" + QString::number(GetTotalCudaMemory() / MB, 'f', 2) + " MB", "memory");

	QObject::connect(&gTransferFunction, SIGNAL(FunctionChanged()), this, SLOT(OnUpdateTransferFunction()));
	QObject::connect(&gTransferFunction, SIGNAL(SettingsChanged()), this, SLOT(OnUpdateTransferFunctionSettings()));
	QObject::connect(&gCamera, SIGNAL(Changed()), this, SLOT(OnUpdateCamera()));
	QObject::connect(&gLighting, SIGNAL(Changed()), this, SLOT(OnUpdateLighting()));
	QObject::connect(&gLighting.Background(), SIGNAL(Changed()), this, SLOT(OnUpdateLighting()));

	QObject::connect(&gStatus, SIGNAL(RenderPause(const bool&)), this, SLOT(OnRenderPause(const bool&)));

	// Try to load appearance/lighting/camera presets with the same name as the loaded file
	gStatus.SetLoadPreset(QFileInfo(m_FileName).baseName());

	// Keep track of frames/second
	CTiming FPS, RenderImage, BlurImage, PostProcessImage, DenoiseImage;

	ResetRenderCanvasView();

	float FpsTime;

	try
	{
		while (!m_Abort)
		{
			if (m_Pause)
				continue;

			if (m_startTesting) {
				Log("Resetting render progress for testing", "control");

				FPS.Reset();
				RenderImage.Reset();
				BlurImage.Reset();
				PostProcessImage.Reset();
				DenoiseImage.Reset();

				gScene.SetNoIterations(0);
				ResetRenderCanvasView();

				m_doTests = true;
				m_startTesting = false;
			}

			gStatus.SetPreRenderFrame();

			// CUDA time for profiling
 			CCudaTimer TmrFps;

			SceneCopy = gScene;

			gStatus.SetStatisticChanged("Camera", "Position", FormatVector(SceneCopy.m_Camera.m_From));
			gStatus.SetStatisticChanged("Camera", "Target", FormatVector(SceneCopy.m_Camera.m_Target));
			gStatus.SetStatisticChanged("Camera", "Up Vector", FormatVector(SceneCopy.m_Camera.m_Up));

			// Resizing the image canvas requires special attention
			if (SceneCopy.m_DirtyFlags.HasFlag(FilmResolutionDirty))
			{
				// Allocate host image buffer, this thread will blit it's frames to this buffer
				free(m_pRenderImage);
				m_pRenderImage = NULL;

				m_pRenderImage = (CColorRgbLdr*)malloc(SceneCopy.m_Camera.m_Film.m_Resolution.GetNoElements() * sizeof(CColorRgbLdr));

				if (m_pRenderImage)
					memset(m_pRenderImage, 0, SceneCopy.m_Camera.m_Film.m_Resolution.GetNoElements() * sizeof(CColorRgbLdr));
			
				gStatus.SetStatisticChanged("Host Memory", "LDR Frame Buffer", QString::number(3 * SceneCopy.m_Camera.m_Film.m_Resolution.GetNoElements() * sizeof(CColorRgbLdr) / MB, 'f', 2), "MB");

				SceneCopy.SetNoIterations(0);


				///temp
				//if (pDevStates != nullptr)
					//HandleCudaError(cudaFree(pDevStates));

				//pDevStates = InitStates(SceneCopy.m_Camera.m_Film.m_Resolution.GetNoElements());
				///temp

				Log("Render canvas resized to: " + QString::number(SceneCopy.m_Camera.m_Film.m_Resolution.GetResX()) + " x " + QString::number(SceneCopy.m_Camera.m_Film.m_Resolution.GetResY()) + " pixels", "application-resize");
			}

			// Restart the rendering when when the camera, lights and render params are dirty
			if (SceneCopy.m_DirtyFlags.HasFlag(CameraDirty | LightsDirty | RenderParamsDirty | TransferFunctionDirty))
			{
				ResetRenderCanvasView();

				// Reset no. iterations
				gScene.SetNoIterations(0);
			}

			// At this point, all dirty flags should have been taken care of, since the flags in the original scene are now cleared
			gScene.m_DirtyFlags.ClearAllFlags();

			SceneCopy.m_DenoiseParams.SetWindowRadius(3.0f);
			SceneCopy.m_DenoiseParams.m_LerpC = 0.33f * (max((float)gScene.GetNoIterations(), 1.0f) * 0.035f);//1.0f - powf(1.0f / (float)gScene.GetNoIterations(), 15.0f);//1.0f - expf(-0.01f * (float)gScene.GetNoIterations());
//			SceneCopy.m_DenoiseParams.m_Enabled = false;

			SceneCopy.m_Camera.Update();

			BindConstants(&SceneCopy);

			BindTransferFunctionOpacity(SceneCopy.m_TransferFunctions.m_Opacity);
			BindTransferFunctionDiffuse(SceneCopy.m_TransferFunctions.m_Diffuse);
			BindTransferFunctionSpecular(SceneCopy.m_TransferFunctions.m_Specular);
			BindTransferFunctionRoughness(SceneCopy.m_TransferFunctions.m_Roughness);
			BindTransferFunctionEmission(SceneCopy.m_TransferFunctions.m_Emission);

			BindRenderCanvasView(SceneCopy.m_Camera.m_Film.m_Resolution);

			//Render(0, SceneCopy, RenderImage, BlurImage, PostProcessImage, DenoiseImage);
  			//Render(SceneCopy.m_AlgorithmType, SceneCopy, RenderImage, BlurImage, PostProcessImage, DenoiseImage);
			Render(SceneCopy, RenderImage, BlurImage, PostProcessImage, DenoiseImage, pDevStates);
		
			gScene.SetNoIterations(gScene.GetNoIterations() + 1);

			gStatus.SetStatisticChanged("Timings", "Render Image", QString::number(RenderImage.m_FilteredDuration, 'f', 2), "ms.");
			gStatus.SetStatisticChanged("Timings", "Blur Estimate", QString::number(BlurImage.m_FilteredDuration, 'f', 2), "ms.");
			gStatus.SetStatisticChanged("Timings", "Post Process Estimate", QString::number(PostProcessImage.m_FilteredDuration, 'f', 2), "ms.");
			gStatus.SetStatisticChanged("Timings", "De-noise Image", QString::number(DenoiseImage.m_FilteredDuration, 'f', 2), "ms.");

			FPS.AddDuration(1000.0f / TmrFps.ElapsedTime());

 			gStatus.SetStatisticChanged("Performance", "FPS", QString::number(FPS.m_FilteredDuration, 'f', 2), "Frames/Sec.");
 			gStatus.SetStatisticChanged("Performance", "No. Iterations", QString::number(SceneCopy.GetNoIterations()), "Iterations");

			HandleCudaError(cudaMemcpy(m_pRenderImage, GetDisplayEstimate(), SceneCopy.m_Camera.m_Film.m_Resolution.GetNoElements() * sizeof(CColorRgbLdr), cudaMemcpyDeviceToHost));
			//some random test
			//HandleCudaError(cudaMemcpy(m_pRenderImage, GetFrameDisplayEstimate(), SceneCopy.m_Camera.m_Film.m_Resolution.GetNoElements() * sizeof(CColorRgbLdr), cudaMemcpyDeviceToHost));

			// total colour intensity
			/*int r = 0, g = 0, b = 0;
			for (int x = 0; x < SceneCopy.m_Camera.m_Film.m_Resolution.GetResX(); x++) {
				for (int y = 0; y < SceneCopy.m_Camera.m_Film.m_Resolution.GetResY(); y++) {
					r += m_pRenderImage[x + x * y].r;
					g += m_pRenderImage[x + x * y].g;
					b += m_pRenderImage[x + x * y].b;
				}
			}
			std::cout << "Intensity: " << r << ", " << g << ", " << b << std::endl;*/

			gFrameBuffer.Set((unsigned char*)m_pRenderImage, SceneCopy.m_Camera.m_Film.GetWidth(), SceneCopy.m_Camera.m_Film.GetHeight());

			if (m_doTests) {
				int index = m_SaveFrames.indexOf(SceneCopy.GetNoIterations());
				if (index >= 0)
				{
					//const QString ImageFilePath = QApplication::applicationDirPath() + "/Output/" + m_SaveBaseName + "_" + QString::number(SceneCopy.GetNoIterations()) + ".png";
					const QString ImageFilePath = m_TestDir + "/images/" + QString::number(SceneCopy.GetNoIterations()) + ".png";

					//Log("Saving test to: " + ImageFilePath, "conrol");
					SaveImage((unsigned char*)m_pRenderImage, SceneCopy.m_Camera.m_Film.m_Resolution.GetResX(), SceneCopy.m_Camera.m_Film.m_Resolution.GetResY(), ImageFilePath);

					m_SaveFrames.removeAt(index);

					QFile measurements(m_TestDir + "/measurements.csv");
					if (measurements.open(QIODevice::WriteOnly | QIODevice::Append)) {
						size_t freeMemory = -1;
						size_t totalMemory = -1;
						HandleCudaError(cudaMemGetInfo(&freeMemory, &totalMemory));
						QTextStream stream(&measurements);
						QStringList data = {
							QString::number(SceneCopy.GetNoIterations()),
							QString::number(RenderImage.m_TotalTime),
							QString::number(BlurImage.m_TotalTime),
							QString::number(PostProcessImage.m_TotalTime),
							QString::number(DenoiseImage.m_TotalTime),
							QString::number((totalMemory - freeMemory) / MB)
						};
						stream << "\r\n" << data.join(",");
						measurements.close();
					}

					if (m_SaveFrames.length() == 0) {
						m_doTests = false;
					}
				}
			}

 			gStatus.SetPostRenderFrame();
		}
	}
	catch (QString* pMessage)
	{
//		Log(*pMessage + ", rendering will be aborted");

		free(m_pRenderImage);
		m_pRenderImage = NULL;

		gStatus.SetRenderEnd();

		return;
	}

	free(m_pRenderImage);
	m_pRenderImage = NULL;

	UnbindDensityBuffer();

	UnbindTransferFunctionOpacity();
	UnbindTransferFunctionDiffuse();
	UnbindTransferFunctionSpecular();
	UnbindTransferFunctionRoughness();
	UnbindTransferFunctionEmission();
	
	FreeRenderCanvasView();

	// Let others know that we have stopped rendering
	gStatus.SetRenderEnd();

	Log("Device memory: " + QString::number(GetUsedCudaMemory() / MB, 'f', 2) + "/" + QString::number(GetTotalCudaMemory() / MB, 'f', 2) + " MB", "memory");

	// Clear the histogram
	gHistogram.Reset();

	ResetDevice();
}

struct DataSpaceRawCoordinate {
	Vec3f coord;
	DataSpaceRawCoordinate(int x, int y, int z) :
		coord(x, y, z)
	{}
};

struct DataSpaceRawSpacing {
	Vec3f coord;
	explicit DataSpaceRawSpacing(Vec3f coord) :
		coord(coord)
	{}
};

struct DataSpaceRawSpaced {
	Vec3f coord;
	DataSpaceRawSpaced(DataSpaceRawCoordinate coord, DataSpaceRawSpacing spacing) :
		coord(coord.coord * spacing.coord)
	{}
};

bool QRenderThread::Load(QString& FileName)
{
	m_FileName = FileName;

	// Create meta image reader
	vtkSmartPointer<vtkMetaImageReader> MetaImageReader = vtkMetaImageReader::New();

	QFileInfo FileInfo(FileName);

	if (!FileInfo.exists())
	{
		Log(QString(QFileInfo(FileName).filePath().replace("//", "/")).toLatin1() + "  does not exist!", QLogger::Critical);
		return false;
	}

	Log(QString("Loading " + QFileInfo(FileName).fileName()).toLatin1());

	// Exit if the reader can't read the file
	if (!MetaImageReader->CanReadFile(m_FileName.toLatin1()))
	{
		Log(QString("Meta image reader can't read file " + QFileInfo(FileName).fileName()).toLatin1(), QLogger::Critical);
		return false;
	}

	MetaImageReader->SetFileName(m_FileName.toLatin1());

	MetaImageReader->Update();

	if (MetaImageReader->GetErrorCode() != vtkErrorCode::NoError)
	{
		Log("Error loading file " + QString(vtkErrorCode::GetStringFromErrorCode(MetaImageReader->GetErrorCode())));
		return false;
	}

	vtkSmartPointer<vtkImageCast> ImageCast = vtkImageCast::New();
	
	Log("Casting volume data type to short", "grid");

	ImageCast->SetInputConnection(MetaImageReader->GetOutputPort());
	ImageCast->SetOutputScalarTypeToShort();
	ImageCast->Update();

	if (ImageCast->GetErrorCode() != vtkErrorCode::NoError)
	{
		Log("vtkImageCast error: " + QString(vtkErrorCode::GetStringFromErrorCode(MetaImageReader->GetErrorCode())));
		return false;
	}
	
	// Volume resolution
	int* pVolumeResolution = ImageCast->GetOutput()->GetExtent();
	gScene.m_Resolution.SetResXYZ(Vec3i(pVolumeResolution[1] + 1, pVolumeResolution[3] + 1, pVolumeResolution[5] + 1));

	Log("Resolution: " + FormatSize(gScene.m_Resolution.GetResXYZ()) + "", "grid");

	DataSpaceRawCoordinate dataResolution(gScene.m_Resolution.GetResX(), gScene.m_Resolution.GetResY(), gScene.m_Resolution.GetResZ());

	// Intensity range
	double* pIntensityRange = ImageCast->GetOutput()->GetScalarRange();
	gScene.m_IntensityRange.SetMin((float)pIntensityRange[0]);
	gScene.m_IntensityRange.SetMax((float)pIntensityRange[1]);

	Log("Intensity range: [" + QString::number(gScene.m_IntensityRange.GetMin()) + ", " + QString::number(gScene.m_IntensityRange.GetMax()) + "]", "grid");

	// Spacing
	double* pSpacing = ImageCast->GetOutput()->GetSpacing();

	gScene.m_Spacing.x = (float)pSpacing[0];
	gScene.m_Spacing.y = (float)pSpacing[1];
	gScene.m_Spacing.z = (float)pSpacing[2];

	Log("Spacing: " + FormatSize(gScene.m_Spacing, 2), "grid");

	DataSpaceRawSpacing dataSpacing(gScene.m_Spacing);

	DataSpaceRawSpaced dataPhysicalSize(dataResolution, dataSpacing);

	// Compute physical size
	const Vec3f PhysicalSize(Vec3f(gScene.m_Spacing.x * (float)gScene.m_Resolution.GetResX(), gScene.m_Spacing.y * (float)gScene.m_Resolution.GetResY(), gScene.m_Spacing.z * (float)gScene.m_Resolution.GetResZ()));

	// Compute the volume's bounding box
	gScene.m_BoundingBox.m_MinP	= Vec3f(0.0f);
	gScene.m_BoundingBox.m_MaxP	= PhysicalSize / PhysicalSize.Max();

	// Define the size of the respective voxel axis in world coordinates
	// +2 is added to the resolution to account for the two extra virtual voxel we add around the existing volume
	// -1 is substracted, because we need to see voxels as fenceposts instead of fences (fencepost problem)
	gScene.m_VoxelSizeWorld = make_float3(
		gScene.m_BoundingBox.m_MaxP.x / (gScene.m_Resolution.GetResX() + 2 - 1),
		gScene.m_BoundingBox.m_MaxP.y / (gScene.m_Resolution.GetResY() + 2 - 1),
		gScene.m_BoundingBox.m_MaxP.z / (gScene.m_Resolution.GetResZ() + 2 - 1));
	
	/*gScene.m_VoxelSizeWorld = make_float3(
		gScene.m_BoundingBox.m_MaxP.x / (gScene.m_Resolution.GetResX() - 1),
		gScene.m_BoundingBox.m_MaxP.y / (gScene.m_Resolution.GetResY() - 1),
		gScene.m_BoundingBox.m_MaxP.z / (gScene.m_Resolution.GetResZ() - 1));*/
	
	Log("Bounding box: " + FormatVector(gScene.m_BoundingBox.m_MinP, 2) + " - " + FormatVector(gScene.m_BoundingBox.m_MaxP), "grid");
	
	const int DensityBufferSize = gScene.m_Resolution.GetNoElements() * sizeof(short);

 	m_pDensityBuffer = (short*)malloc(DensityBufferSize);
  	memcpy(m_pDensityBuffer, ImageCast->GetOutput()->GetScalarPointer(), DensityBufferSize);

	// Gradient magnitude volume
	vtkSmartPointer<vtkImageGradientMagnitude> GradientMagnitude = vtkImageGradientMagnitude::New();
	
	Log("Creating gradient magnitude volume", "grid");
		
	GradientMagnitude->SetDimensionality(3);
	GradientMagnitude->SetInputConnection(ImageCast->GetOutputPort());
	GradientMagnitude->Update();

	vtkImageData* GradientMagnitudeBuffer = GradientMagnitude->GetOutput();
	
	// Scalar range of the gradient magnitude
	double* pGradientMagnitudeRange = GradientMagnitudeBuffer->GetScalarRange();
	
	gScene.m_GradientMagnitudeRange.SetMin((float)pGradientMagnitudeRange[0]);
	gScene.m_GradientMagnitudeRange.SetMax((float)pGradientMagnitudeRange[1]);
	
	Log("Gradient magnitude range: [" + QString::number(gScene.m_GradientMagnitudeRange.GetMin(), 'f', 2) + " - " + QString::number(gScene.m_GradientMagnitudeRange.GetMax(), 'f', 2) + "]", "grid");
	
	const int GradientMagnitudeBufferSize = gScene.m_Resolution.GetNoElements() * sizeof(short);
	
	m_pGradientMagnitudeBuffer = (short*)malloc(GradientMagnitudeBufferSize);
	memcpy(m_pGradientMagnitudeBuffer, GradientMagnitudeBuffer->GetScalarPointer(), GradientMagnitudeBufferSize);

	// Build the histogram
	Log("Creating gradient magnitude histogram", "grid");

	vtkSmartPointer<vtkImageAccumulate> GradMagHistogram = vtkSmartPointer<vtkImageAccumulate>::New();

	GradMagHistogram->SetInputConnection(GradientMagnitude->GetOutputPort());
	GradMagHistogram->SetComponentExtent(0, 255, 0, 0, 0, 0);
	GradMagHistogram->SetComponentOrigin(0, 0, 0);
	GradMagHistogram->SetComponentSpacing(gScene.m_GradientMagnitudeRange.GetRange() / 256.0f, 0, 0);
//	GradMagHistogram->IgnoreZeroOn();
	GradMagHistogram->Update();

	gScene.m_GradMagMean = (float)GradMagHistogram->GetMean()[0];
	gScene.m_GradientFactor = gScene.m_GradMagMean;

	Log("Mean gradient magnitude: " + QString::number(gScene.m_GradMagMean, 'f', 2), "grid");

	Log("Creating density histogram", "grid");

	// Build the histogram
	vtkSmartPointer<vtkImageAccumulate> Histogram = vtkSmartPointer<vtkImageAccumulate>::New();

	Log("Creating histogram", "grid");
 	Histogram->SetInputConnection(ImageCast->GetOutputPort());
 	Histogram->SetComponentExtent(0, 255, 0, 0, 0, 0);
 	Histogram->SetComponentOrigin(gScene.m_IntensityRange.GetMin(), 0, 0);
 	Histogram->SetComponentSpacing((gScene.m_IntensityRange.GetRange() + 0.01) / 256.0f, 0, 0);
// 	Histogram->IgnoreZeroOn();
 	Histogram->Update();

	// Update the histogram in the transfer function
	gHistogram.SetBins((vtkIdType*)Histogram->GetOutput()->GetScalarPointer(), 256);
	
	gStatus.SetStatisticChanged("Volume", "File", QFileInfo(m_FileName).fileName(), "");
	gStatus.SetStatisticChanged("Volume", "Bounding Box", "", "");
	gStatus.SetStatisticChanged("Bounding Box", "Min", FormatVector(gScene.m_BoundingBox.m_MinP, 2), "m");
	gStatus.SetStatisticChanged("Bounding Box", "Max", FormatVector(gScene.m_BoundingBox.m_MaxP, 2), "m");
	gStatus.SetStatisticChanged("Volume", "Physical Size", FormatSize(PhysicalSize, 2), "mm");
	gStatus.SetStatisticChanged("Volume", "Resolution", FormatSize(gScene.m_Resolution.GetResXYZ()), "Voxels");
	gStatus.SetStatisticChanged("Volume", "Spacing", FormatSize(gScene.m_Spacing, 2), "mm");
	gStatus.SetStatisticChanged("Volume", "No. Voxels", QString::number(gScene.m_Resolution.GetNoElements()), "Voxels");
	gStatus.SetStatisticChanged("Volume", "Density Range", "[" + QString::number(gScene.m_IntensityRange.GetMin()) + ", " + QString::number(gScene.m_IntensityRange.GetMax()) + "]", "");



	// Try to save our own file
	// adapt path !
	std::string filePath = "../exposure-render.release110/Source/Examples/tunnelAligned.mhd";
	std::string filePathRaw = "../exposure-render.release110/Source/Examples/tunnelAligned.raw";

	struct stat buffer;
	if (!((stat(filePath.c_str(), &buffer) == 0))) {
		
		// Create an 3D cross
		/*const int parts = 7;
		const int partSize = 2; //change this value for bigger volumes
		const int width = parts * partSize;
		const int height = parts * partSize;
		const int depth = parts * partSize;

		std::unique_ptr<short[]> img(new short[width * height * depth]);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;
					img[id] = 0;
					if ((row >=  3 * partSize && row < 4 * partSize)
						|| (col >= 3 * partSize && col < 4 * partSize)
						|| (dep >= 3 * partSize && dep < 4 * partSize)
							) {
						img[id] = 100;
					}
					if ((dep >= 3 * partSize && dep < 4 * partSize)
						&& (row >= 5 * partSize && row < 6 * partSize)
						&& ((col >= 1 * partSize && col < 2 * partSize) || (col >= 5 * partSize && col < 6 * partSize))
						) {
						img[id] = 50;
					}
				}
			}
		}*/

		// gradient
		/*const int N = 10;
		const int width = N;
		const int height = 4;
		const int depth = 4;

		std::unique_ptr<short[]> img(new short[width * height * depth]);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < N; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;
					img[id] = col + 1;
					if (row == 0 || row == height - 1
						|| dep == 0 || dep == depth - 1)
						img[id] = 0;
				}
			}
		}*/

		// Tunnel Grid Aligned
		const int tunnelDiameter = 16;
		const int tunnelLengths = 32;
		const int buffer = 16;

		const int width = tunnelLengths + tunnelDiameter * 2 + buffer * 2;
		const int height = width;
		const int depth = 2 * tunnelLengths + tunnelDiameter + buffer;

		std::unique_ptr<short[]> img(new short[width * height * depth]);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;

					img[id] = 0;

					if (dep < buffer) {
						img[id] = 1000;
						if (col >= buffer && col < buffer + tunnelDiameter + (buffer - dep)) {
							img[id] = 100;
							if (row < buffer)
								img[id] = 600;
						}
					}
					else if (dep < buffer + tunnelLengths) {
						img[id] = 1000;
						if (col >= buffer && col < buffer + tunnelDiameter) {
							img[id] = 100;
							if (row < buffer)
								img[id] = 600;
						}
					}
					else if (dep < buffer + tunnelLengths + tunnelDiameter) {
						img[id] = 1000;
						if (col >= buffer && col < buffer + tunnelLengths + 2 * tunnelDiameter) {
							img[id] = 100;
							if (row < buffer)
								img[id] = 600;
						}
					}
					else {
						img[id] = 1000;
						if (col >= buffer + tunnelLengths + tunnelDiameter && col < buffer + tunnelLengths + 2 * tunnelDiameter) {
							img[id] = 100;
							if (row < buffer)
								img[id] = 600;
						}
					}
				}
			}
		}
		

		// Create an Cornell box image
		/*const int size = 128;
		const int width = size;
		const int height = size;
		const int depth = size;
		int wt = 8;

		short* img;
		img = (short*)malloc(width*height*depth*sizeof(short));
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;
					
					img[id] = 0;

					if (dep > 0) {
						const bool wall = (row < wt || row >= size - wt) || (col < wt || col >= size - wt) || dep >= size - wt;

						if (wall)
							img[id] = 100;

						const bool floor = row < wt;

						if (floor)
							img[id] = 50;
					}
				}
			}
		}*/

		/*
		// Intensity gradient
		const int size = 128;
		const int width = size;
		const int height = size;
		const int depth = size;

		int currentSum = 0;
		int density = -5;
		int itt = 0;

		int target = pow(2, itt);

		std::unique_ptr<short[]> img(new short[width * height * depth]);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;

					img[id] = density;
					currentSum++;

					if (currentSum == target) {
						itt++;
						currentSum = 0;
						density++;
						target = pow(2, itt);
					}
				}
			}
		}*/
		/*
		//Tunnel
		const int size = 64;
		const int width = size;
		const int height = size;

		const int tunnelRadius = size / 4;
		const int wallThickness = size / 8;

		const int depth = size + wallThickness;

		std::unique_ptr<short[]> img(new short[width * height * depth]);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;

					img[id] = 0;

					if (dep < wallThickness) {
						img[id] = 1000;
						if (col >= wallThickness && col < wallThickness + tunnelRadius) {
							img[id] = 100;
							if (row < wallThickness)
								img[id] = 1000;
						}
					}
					else {
						int c = dep - 8;
						if (c < (depth - wallThickness) / 2) {
							if (col >= c)
								img[id] = 500;
							if (col >= c + wallThickness)
								img[id] = 100;
							if (col >= c + wallThickness + tunnelRadius)
								img[id] = 500;
							if (col >= c + 2 * wallThickness + tunnelRadius)
								img[id] = 0;
						}
						else {
							c = depth - dep;
							if (col >= c)
								img[id] = 500;
							if (col >= c + wallThickness)
								img[id] = 100;
							if (col >= c + wallThickness + tunnelRadius)
								img[id] = 500;
							if (col >= c + 2 * wallThickness + tunnelRadius)
								img[id] = 0;
						}
						if (row < wallThickness)
							img[id] = 500;
					}
				}
			}
		}
		*/
		
		/*
		//Checkerbox
		const int size = 128;
		const int width = size;
		const int height = size;
		const int depth = size;

		std::unique_ptr<short[]> img(new short[width * height * depth]);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;

					img[id] = ((row+col) % 2) * 1000;
					if (row + col < size)
						img[id] = 0;
				}
			}
		}
		*/

		/*
		//Chekerbox
		const int size = 2;
		const int width = size;
		const int height = size;
		const int depth = size;

		std::unique_ptr<short[]> img(new short[width * height * depth]);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				for (int dep = 0; dep < depth; dep++) {
					int id = col + row * width + dep * width * height;

					img[id] = ((row + col) % 2) * 1000;
				}
			}
		}
		*/

		// Convert the c-style image to a vtkImageData
		vtkSmartPointer<vtkImageImport> imageImport =
			vtkSmartPointer<vtkImageImport>::New();
		imageImport->SetDataSpacing(1, 1, 1);
		imageImport->SetDataOrigin(0, 0, 0);
		imageImport->SetWholeExtent(0, width - 1, 0, height - 1, 0, depth - 1);
		imageImport->SetDataExtentToWholeExtent();
		imageImport->SetDataScalarTypeToShort();
		imageImport->SetNumberOfScalarComponents(1);
		imageImport->SetImportVoidPointer(img.get());
		imageImport->Update();

		vtkSmartPointer<vtkImageCast> castFilter =
			vtkSmartPointer<vtkImageCast>::New();
		castFilter->SetOutputScalarTypeToShort();
		castFilter->SetInputConnection(imageImport->GetOutputPort());
		castFilter->Update();

		vtkSmartPointer<vtkMetaImageWriter> writer =
			vtkSmartPointer<vtkMetaImageWriter>::New();
		writer->SetInputConnection(castFilter->GetOutputPort());
		writer->SetFileName(filePath.c_str());
		writer->SetRAWFileName(filePathRaw.c_str());
		writer->Write();
	}
	
	return true;
}

void QRenderThread::OnUpdateTransferFunction(void)
{
	QMutexLocker Locker(&gSceneMutex);

	QTransferFunction TransferFunction = gTransferFunction;

	gScene.m_TransferFunctions.m_Opacity.m_NoNodes		= TransferFunction.GetNodes().size();
	gScene.m_TransferFunctions.m_Diffuse.m_NoNodes		= TransferFunction.GetNodes().size();
	gScene.m_TransferFunctions.m_Specular.m_NoNodes		= TransferFunction.GetNodes().size();
	gScene.m_TransferFunctions.m_Emission.m_NoNodes = TransferFunction.GetNodes().size();
	gScene.m_TransferFunctions.m_Roughness.m_NoNodes = TransferFunction.GetNodes().size();

	for (int i = 0; i < TransferFunction.GetNodes().size(); i++)
	{
		QNode& Node = TransferFunction.GetNode(i);

		const float Intensity = Node.GetIntensity();

		// Positions
		gScene.m_TransferFunctions.m_Opacity.m_P[i] = Intensity;
		gScene.m_TransferFunctions.m_Diffuse.m_P[i] = Intensity;
		gScene.m_TransferFunctions.m_Specular.m_P[i] = Intensity;
		gScene.m_TransferFunctions.m_Emission.m_P[i] = Intensity;
		gScene.m_TransferFunctions.m_Roughness.m_P[i] = Intensity;

		// Colors
		gScene.m_TransferFunctions.m_Opacity.m_C[i] = CColorRgbHdr(Node.GetOpacity());
		gScene.m_TransferFunctions.m_Diffuse.m_C[i] = CColorRgbHdr(Node.GetDiffuse().redF(), Node.GetDiffuse().greenF(), Node.GetDiffuse().blueF());
		gScene.m_TransferFunctions.m_Specular.m_C[i] = CColorRgbHdr(Node.GetSpecular().redF(), Node.GetSpecular().greenF(), Node.GetSpecular().blueF());
		gScene.m_TransferFunctions.m_Emission.m_C[i] = 500.0f * CColorRgbHdr(Node.GetEmission().redF(), Node.GetEmission().greenF(), Node.GetEmission().blueF());

		const float Roughness = 1.0f - expf(-Node.GetGlossiness());

		gScene.m_TransferFunctions.m_Roughness.m_C[i] = CColorRgbHdr(Roughness * 250.0f);
	}

	gScene.m_DensityScale = TransferFunction.GetDensityScale();
	gScene.m_GradientFactor = TransferFunction.GetGradientFactor();

	gScene.m_DirtyFlags.SetFlag(TransferFunctionDirty);
}

void QRenderThread::OnUpdateTransferFunctionSettings(void) {
	QMutexLocker Locker(&gSceneMutex);

	QTransferFunction TransferFunction = gTransferFunction;

	gScene.m_ShadingType = TransferFunction.GetShadingType();
	gScene.m_AlgorithmType = TransferFunction.GetAlgorithmType();
	gScene.m_ScatterType = TransferFunction.GetScatterType();
	gScene.m_MaxBounces = TransferFunction.GetNrOfBounces();
	gScene.m_PostProcessingSteps = TransferFunction.GetPostProcessingSteps();
	gScene.m_StepSizeFactor = TransferFunction.GetPrimaryStepSize();
	gScene.m_StepSizeFactorShadow = TransferFunction.GetSecondarStepSize();
	gScene.m_ScatteringHeadstart = TransferFunction.GetScatteringHeadstart();

	if ((gScene.m_AlgorithmType == SINGLE_SCATTERING || gScene.m_AlgorithmType == MULTIPLE_SCATTERING)
		&& m_CurFilterMode != cudaFilterModeLinear) {
		cudaExtent Res;
		Res.width = gScene.m_Resolution[0];
		Res.height = gScene.m_Resolution[1];
		Res.depth = gScene.m_Resolution[2];

		// Bind density buffer to texture
		Log("Copying density volume to device with FilterModeLinear", "grid");

		UnbindDensityBuffer();
		BindDensityBuffer((short*)m_pDensityBuffer, Res, cudaFilterModeLinear);
		m_CurFilterMode = cudaFilterModeLinear;
	}

	if (gScene.m_AlgorithmType == 2) {
		InitPreCalculated();
	}
	else if (gScene.m_AlgorithmType == 3) {
		InitOpacityGradient(CScene(gScene));
	}
	else if (gScene.m_AlgorithmType == 4) {
		InitFloodFill();
	}
	else if (gScene.m_AlgorithmType == PROPERTY_BASED && m_CurFilterMode != cudaFilterModePoint) {
		cudaExtent Res;
		Res.width = gScene.m_Resolution[0];
		Res.height = gScene.m_Resolution[1];
		Res.depth = gScene.m_Resolution[2];

		// Bind density buffer to texture
		Log("Copying density volume to device with FilterModePoint", "grid");

		UnbindDensityBuffer();
		BindDensityBuffer((short*)m_pDensityBuffer, Res, cudaFilterModePoint);
		m_CurFilterMode = cudaFilterModePoint;
	}

	gScene.m_DirtyFlags.SetFlag(TransferFunctionDirty);
}

void QRenderThread::InitPreCalculated() {
	CScene SceneCopy = gScene;
	InitPreCalculatedCore(SceneCopy, m_pDensityBuffer);
}

int* QRenderThread::InitFloodFill() {
	
	// using lambda to compare elements.
	auto compare = [](Vec4i lhs, Vec4i rhs)
	{
		return lhs.w > rhs.w;
	};

	int** results;
	//results = (int**)malloc(sizeof(int*) * gScene.m_Lighting.m_NoLights);
	results = (int**)malloc(sizeof(int*) * 1);

	for (int lightIndex = 0; lightIndex < gScene.m_Lighting.m_NoLights; lightIndex++) {
		// Make an array to store our results
		int* result;
		result = (int*)malloc(gScene.m_Resolution.GetNoElements() * sizeof(int));
		for (int i = 0; i < gScene.m_Resolution.GetNoElements(); i++) {
			result[i] = -1;
		}
		//results[lightIndex] = result;
		results[0] = result;

		// Make a new priority queue
		std::priority_queue<Vec4i, std::vector<Vec4i>, decltype(compare)> queue(compare);

		// Setup regarding the current light
		CLight light = gScene.m_Lighting.m_Lights[lightIndex];
		std::cout << "LightIndex: " << lightIndex << ", Type: " << light.m_T << std::endl;
		if (light.m_T == 0) { // Area light
			// Get the direction from the light center point to the middle of the volume
			Vec3f dir = Vec3f(0.5) - light.m_P;
			dir.Normalize();

			float nearest;
			float furthest;
			IntersectBox(light.m_P, dir, &nearest, &furthest);

			Vec3f entryWorld = light.m_P + dir * nearest;
			Vec3f volExit = light.m_P + dir * furthest;

			Vec3f* points = new Vec3f[4];
			points[0] = light.m_P + (-light.m_HalfWidth * light.m_U) + (-light.m_HalfHeight * light.m_V);
			points[1] = light.m_P + (+light.m_HalfWidth * light.m_U) + (-light.m_HalfHeight * light.m_V);
			points[2] = light.m_P + (+light.m_HalfWidth * light.m_U) + (+light.m_HalfHeight * light.m_V);
			points[3] = light.m_P + (-light.m_HalfWidth * light.m_U) + (+light.m_HalfHeight * light.m_V);
			/*
			Vec4i* pointsVol = new Vec4i[4];
			for (int i = 0; i < 4; i++) {
				IntersectBox(points[i], volExit - points[i], &nearest, &furthest);
				Vec3f entry = points[i] + (volExit - points[i]) * nearest;
				pointsVol[i] = Vec4i(
					Clamp(entry.x / gScene.m_VoxelSizeWorld.x, 0, gScene.m_Resolution.GetResX() - 1),
					Clamp(entry.y / gScene.m_VoxelSizeWorld.y, 0, gScene.m_Resolution.GetResY() - 1),
					Clamp(entry.z / gScene.m_VoxelSizeWorld.z, 0, gScene.m_Resolution.GetResZ() - 1),
					0
				);
			}*/
			float threshold = 0.5;
			for (int x = 0; x < gScene.m_Resolution.GetResX(); x++) {
				for (int y = 0; y < gScene.m_Resolution.GetResY(); y++) {
					for (int z = 0; z < gScene.m_Resolution.GetResZ(); z++) {
						// Check if this is a point along the edge
						if (x == 0 || x == gScene.m_Resolution.GetResX() - 1
							|| y == 0 || y == gScene.m_Resolution.GetResY() - 1
							|| z == 0 || z == gScene.m_Resolution.GetResZ() - 1) {

							Vec3f p = Vec3f(x * gScene.m_VoxelSizeWorld.x, y * gScene.m_VoxelSizeWorld.y, z * gScene.m_VoxelSizeWorld.z);
							float maxDot = -1;
							for (int i = 0; i < 4; i++) {
								Vec3f pToL = points[i] - p;
								pToL.Normalize();
								if (x == 0) {
									maxDot = Fmaxf(maxDot, Vec3f(-1, 0, 0).Dot(pToL));
									maxDot = Fmaxf(maxDot, Vec3f(-1, 0, 0).Dot(-light.m_N));
								}
								if (x == gScene.m_Resolution.GetResX() - 1) {
									maxDot = Fmaxf(maxDot, Vec3f(1, 0, 0).Dot(pToL));
									maxDot = Fmaxf(maxDot, Vec3f(1, 0, 0).Dot(-light.m_N));
								}
								if (y == 0) {
									maxDot = Fmaxf(maxDot, Vec3f(0, -1, 0).Dot(pToL));
									maxDot = Fmaxf(maxDot, Vec3f(0, -1, 0).Dot(-light.m_N));
								}
								if (y == gScene.m_Resolution.GetResY() - 1) {
									maxDot = Fmaxf(maxDot, Vec3f(0, 1, 0).Dot(pToL));
									maxDot = Fmaxf(maxDot, Vec3f(0, 1, 0).Dot(-light.m_N));
								}
								if (z == 0) {
									maxDot = Fmaxf(maxDot, Vec3f(0, 0, -1).Dot(pToL));
									maxDot = Fmaxf(maxDot, Vec3f(0, 0, -1).Dot(-light.m_N));
								}
								if (z == gScene.m_Resolution.GetResZ() - 1) {
									maxDot = Fmaxf(maxDot, Vec3f(0, 0, 1).Dot(pToL));
									maxDot = Fmaxf(maxDot, Vec3f(0, 0, 1).Dot(-light.m_N));
								}

							}

							if (maxDot > threshold) {
								Vec4i start = Vec4i(x, y, z, 0);
								int pos1D = start.x + gScene.m_Resolution.GetResX() * start.y + gScene.m_Resolution.GetResX() * gScene.m_Resolution.GetResY() * start.z;
								float normalizedDensity = (m_pDensityBuffer[pos1D] - gScene.m_IntensityRange.GetMin()) / gScene.m_IntensityRange.GetRange();
								float opacity = gScene.m_TransferFunctions.m_Opacity.F(normalizedDensity).r;
								start.w = (int)(opacity * 100 + 1) * (2 - (maxDot + 1));
								result[pos1D] = start.w;
								queue.push(start);
							}
						}
					}
				}
			}
		}
		else if (light.m_T == 1) { // Background light
			// lets start with background only for now
			continue;

			for (int x = 0; x < gScene.m_Resolution.GetResX(); x++) {
				for (int y = 0; y < gScene.m_Resolution.GetResY(); y++) {
					for (int z = 0; z < gScene.m_Resolution.GetResZ(); z++) {
						// Check if this is a point along the edge
						if (x == 0 || x == gScene.m_Resolution.GetResX() - 1
							|| y == 0 || y == gScene.m_Resolution.GetResY() - 1
							|| z == 0 || z == gScene.m_Resolution.GetResZ() - 1) {
							Vec4i start = Vec4i(x, y, z, 0);
							int pos1D = start.x + gScene.m_Resolution.GetResX() * start.y + gScene.m_Resolution.GetResX() * gScene.m_Resolution.GetResY() * start.z;
							float normalizedDensity = (m_pDensityBuffer[pos1D] - gScene.m_IntensityRange.GetMin()) / gScene.m_IntensityRange.GetRange();
							float opacity = gScene.m_TransferFunctions.m_Opacity.F(normalizedDensity).r;
							start.w = (int)(opacity * 100 + 1);
							result[pos1D] = start.w;
							queue.push(start);
						}
					}
				}
			}
		}
		else { // Unknown type. skip
			continue;
		}

		int i = 1;

		while (!queue.empty()) {
			Vec4i point = queue.top();
			queue.pop();
			//std::cout << "Point: " << point.x << ", " << point.y << ", " << point.z << ", " << point.w << ". Result: " << result[pos1D] << std::endl;

			if (i++ % (gScene.m_Resolution.GetNoElements() / 100) == 0)
				std::cout << "i: " << i << " of " << gScene.m_Resolution.GetNoElements()
				<< " = " << ((int)((float)i / gScene.m_Resolution.GetNoElements() * 10000) / 100.0f) 
				<< "%. Point: " << point.x << ", " << point.y << ", " << point.z << ". Value: " << point.w << std::endl;

			for (int x = -1; x < 2; x++) {
				for (int y = -1; y < 2; y++) {
					for (int z = -1; z < 2; z++) {
						Vec4i newPoint = Vec4i(point.x + x, point.y + y, point.z + z, point.w);
						int pos1D = newPoint.x + gScene.m_Resolution.GetResX() * newPoint.y + gScene.m_Resolution.GetResX() * gScene.m_Resolution.GetResY() * newPoint.z;
						//std::cout << "newPoint: " << point.x << ", " << point.y << ", " << point.z << ", " << point.w << ". Result: " << result[pos1D] << std::endl;
						// Check if we have a usable point. Otherwise continue to the next point
						if (newPoint.x < 0 || newPoint.x >= gScene.m_Resolution.GetResX()
							|| newPoint.y < 0 || newPoint.y >= gScene.m_Resolution.GetResY()
							|| newPoint.z < 0 || newPoint.z >= gScene.m_Resolution.GetResZ()
							|| result[pos1D] != -1)
							continue;

						float normalizedDensity = (m_pDensityBuffer[pos1D] - gScene.m_IntensityRange.GetMin()) / gScene.m_IntensityRange.GetRange();
						float opacity = gScene.m_TransferFunctions.m_Opacity.F(normalizedDensity).r;

						// Some function to map opacity to step size. Could emulating stepsizes of woodcock more closely help?
						newPoint.w += (int)(opacity * 100 + 1);

						// Store the result in an array
						result[pos1D] = newPoint.w;

						// Add new point to the queue
						queue.push(newPoint);
					}
				}
			}
		}
	}

	/*for (int i = 0; i < gScene.m_Resolution.GetNoElements(); i++) {
		//results[0][i] = 999;
		std::cout << "index: " << i << ", result: " + std::to_string(results[0][i]) << std::endl;
	}*/

	return results[0];
}

/// <summary>
/// Checks of the ray intersects with the boundingbox of the volume and fills minT and maxT with the correct values
/// </summary>
bool QRenderThread::IntersectBox(const Vec3f Pe, const Vec3f& Dir, float* pNearT, float* pFarT)
{
	const Vec3f InvR = Vec3f(1.0f, 1.0f, 1.0f) / Dir;
	const Vec3f BottomT = InvR * (gScene.m_BoundingBox.GetMinP() - Pe);
	const Vec3f TopT = InvR * (gScene.m_BoundingBox.GetMaxP() - Pe);
	const Vec3f MinT = MinVec3f(TopT, BottomT);
	const Vec3f MaxT = MaxVec3f(TopT, BottomT);
	const float LargestMinT = fmaxf(fmaxf(MinT.x, MinT.y), fmaxf(MinT.x, MinT.z));
	const float LargestMaxT = fminf(fminf(MaxT.x, MaxT.y), fminf(MaxT.x, MaxT.z));

	*pNearT = LargestMinT;
	*pFarT = LargestMaxT;

	return LargestMaxT > LargestMinT;
}

void QRenderThread::OnUpdateCamera(void)
{
	QMutexLocker Locker(&gSceneMutex);

	gScene.m_Camera.m_Film.m_Exposure = 1.0f - gCamera.GetFilm().GetExposure();

	if (gCamera.GetFilm().IsDirty())
	{
		const int FilmWidth	= gCamera.GetFilm().GetWidth();
		const int FilmHeight = gCamera.GetFilm().GetHeight();

 		gScene.m_Camera.m_Film.m_Resolution.SetResX(FilmWidth);
		gScene.m_Camera.m_Film.m_Resolution.SetResY(FilmHeight);
		gScene.m_Camera.Update();
		gCamera.GetFilm().UnDirty();
// 		// 
 		gScene.m_DirtyFlags.SetFlag(FilmResolutionDirty);
	}

// 	gScene.m_Camera.m_From	= gCamera.GetFrom();
// 	gScene.m_Camera.m_Target	= gCamera.GetTarget();
// 	gScene.m_Camera.m_Up		= gCamera.GetUp();

	gScene.m_Camera.Update();

	// Aperture
	gScene.m_Camera.m_Aperture.m_Size	= gCamera.GetAperture().GetSize();

	// Projection
	gScene.m_Camera.m_FovV = gCamera.GetProjection().GetFieldOfView();

	// Focus
	gScene.m_Camera.m_Focus.m_Type			= (CFocus::EType)gCamera.GetFocus().GetType();
	gScene.m_Camera.m_Focus.m_FocalDistance = gCamera.GetFocus().GetFocalDistance();

	gScene.m_DenoiseParams.m_Enabled = gCamera.GetFilm().GetNoiseReduction();

	gScene.m_DirtyFlags.SetFlag(CameraDirty);
}

void QRenderThread::OnUpdateLighting(void)
{
	QMutexLocker Locker(&gSceneMutex);

	gScene.m_Lighting.Reset();

	if (gLighting.Background().GetEnabled())
	{
		CLight BackgroundLight;

		BackgroundLight.m_T	= 1;

		BackgroundLight.m_ColorTop		= gLighting.Background().GetIntensity() * CColorRgbHdr(gLighting.Background().GetTopColor().redF(), gLighting.Background().GetTopColor().greenF(), gLighting.Background().GetTopColor().blueF());
		BackgroundLight.m_ColorMiddle	= gLighting.Background().GetIntensity() * CColorRgbHdr(gLighting.Background().GetMiddleColor().redF(), gLighting.Background().GetMiddleColor().greenF(), gLighting.Background().GetMiddleColor().blueF());
		BackgroundLight.m_ColorBottom	= gLighting.Background().GetIntensity() * CColorRgbHdr(gLighting.Background().GetBottomColor().redF(), gLighting.Background().GetBottomColor().greenF(), gLighting.Background().GetBottomColor().blueF());
		
		BackgroundLight.Update(gScene.m_BoundingBox);

		gScene.m_Lighting.AddLight(BackgroundLight);

		// TODO: remove printing of text
		/*
		std::cout << "name: backgroundTop, Intensity: " << gLighting.Background().GetIntensity() 
			<< ", Color(RGB): " << gLighting.Background().GetTopColor().redF() << ", " << gLighting.Background().GetTopColor().greenF() << ", " << gLighting.Background().GetTopColor().blueF() 
			<< ", ColorFINAL(RGB): " << BackgroundLight.m_ColorTop.r << ", " << BackgroundLight.m_ColorTop.g << ", " << BackgroundLight.m_ColorTop.b << std::endl;

		std::cout << "name: backgroundMiddle, Intensity: " << gLighting.Background().GetIntensity()
			<< ", Color(RGB): " << gLighting.Background().GetMiddleColor().redF() << ", " << gLighting.Background().GetMiddleColor().greenF() << ", " << gLighting.Background().GetMiddleColor().blueF()
			<< ", ColorFINAL(RGB): " << BackgroundLight.m_ColorMiddle.r << ", " << BackgroundLight.m_ColorMiddle.g << ", " << BackgroundLight.m_ColorMiddle.b << std::endl;

		std::cout << "name: backgroundBottom, Intensity: " << gLighting.Background().GetIntensity()
			<< ", Color(RGB): " << gLighting.Background().GetBottomColor().redF() << ", " << gLighting.Background().GetBottomColor().greenF() << ", " << gLighting.Background().GetBottomColor().blueF()
			<< ", ColorFINAL(RGB): " << BackgroundLight.m_ColorBottom.r << ", " << BackgroundLight.m_ColorBottom.g << ", " << BackgroundLight.m_ColorBottom.b << std::endl;
			*/
	}

	for (int i = 0; i < gLighting.GetLights().size(); i++)
	{
		QLight& Light = gLighting.GetLights()[i];

		CLight AreaLight;

		AreaLight.m_T			= 0;
		AreaLight.m_Theta		= Light.GetTheta() / RAD_F;
		AreaLight.m_Phi			= Light.GetPhi() / RAD_F;
		AreaLight.m_Width		= Light.GetWidth();
		AreaLight.m_Height		= Light.GetHeight();
		AreaLight.m_Distance	= Light.GetDistance();
		AreaLight.m_Color		= Light.GetIntensity() * CColorRgbHdr(Light.GetColor().redF(), Light.GetColor().greenF(), Light.GetColor().blueF());

		//std::cout << "name: " << Light.GetName().toUtf8().constData() << ", Intensity: " << Light.GetIntensity() << ", Color(RGB): " << Light.GetColor().redF() << ", " << Light.GetColor().greenF() << ", " << Light.GetColor().blueF() << ", ColorFINAL(RGB): " << AreaLight.m_Color.r << ", " << AreaLight.m_Color.g << ", " << AreaLight.m_Color.b << std::endl;

		AreaLight.Update(gScene.m_BoundingBox);

		gScene.m_Lighting.AddLight(AreaLight);
	}

	if (gScene.m_ScatterType == 3) {
		cudaExtent ext;
		ext.width = gScene.m_Resolution[0];
		ext.height = gScene.m_Resolution[1];
		ext.depth = gScene.m_Resolution[2];

		int* res = InitFloodFill();
		/*for (int i = 0; i < gScene.m_Resolution.GetNoElements(); i++) {
		std::cout << "index: " << i << ", result: " + std::to_string(res[i]) << std::endl;
		}*/
		UnbindLightPathsBuffer();
		BindLightPathsBuffer(res, ext);
		//free(res);
	}

	gScene.m_DirtyFlags.SetFlag(LightsDirty);
}

void QRenderThread::OnRenderPause(const bool& Pause)
{
	m_Pause = Pause;
}

CColorRgbLdr* QRenderThread::GetRenderImage(void) const
{
	return m_pRenderImage;
}

void StartRenderThread(QString& FileName)
{
	// Create new render thread
 	gpRenderThread = new QRenderThread(FileName);

	// Load the volume
 	if (!gpRenderThread->Load(FileName))
 		return;
 
	// Start the render thread
	gpRenderThread->start();
}

void KillRenderThread(void)
{
 	if (!gpRenderThread)
 		return;
 
	// Kill the render thread
	gpRenderThread->Close();

	// Wait for thread to end
	gpRenderThread->wait();

	// Remove the render thread
	delete gpRenderThread;
	gpRenderThread = NULL;
}

void QRenderThread::CreateIlluminanceTexture() {
	CScene SceneCopy = gScene;

	int SizeIlluminance = SceneCopy.m_Resolution.GetNoElements() * sizeof(float);
	float* pIlluminanceTexture;
	pIlluminanceTexture = (float *)malloc(SizeIlluminance);

	CreateIlluminanceTextureCore(SceneCopy, pIlluminanceTexture);

	// Store Illuminance texture
	unsigned __int64 now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	std::string filePath = "../exposure-render.release110/Source/Examples/illuminance/volume" + std::to_string(now) + ".mhd";
	std::string filePathRaw = "../exposure-render.release110/Source/Examples/illuminance/volume" + std::to_string(now) + ".raw";
	std::chrono::system_clock::now().time_since_epoch();
	
	// Convert the c-style image to a vtkImageData
	vtkSmartPointer<vtkImageImport> imageImport =
		vtkSmartPointer<vtkImageImport>::New();
	imageImport->SetDataSpacing(1, 1, 1);
	imageImport->SetDataOrigin(0, 0, 0);
	imageImport->SetWholeExtent(0, SceneCopy.m_Resolution.GetResX() - 1, 0, SceneCopy.m_Resolution.GetResY() - 1, 0, SceneCopy.m_Resolution.GetResZ() - 1);
	imageImport->SetDataExtentToWholeExtent();
	imageImport->SetDataScalarTypeToFloat();
	imageImport->SetNumberOfScalarComponents(1);
	imageImport->SetImportVoidPointer(pIlluminanceTexture);
	imageImport->Update();

	vtkSmartPointer<vtkImageCast> castFilter =
		vtkSmartPointer<vtkImageCast>::New();
	castFilter->SetOutputScalarTypeToFloat();
	castFilter->SetInputConnection(imageImport->GetOutputPort());
	castFilter->Update();

	vtkSmartPointer<vtkMetaImageWriter> writer =
		vtkSmartPointer<vtkMetaImageWriter>::New();
	writer->SetInputConnection(castFilter->GetOutputPort());
	writer->SetFileName(filePath.c_str());
	writer->SetRAWFileName(filePathRaw.c_str());
	writer->Write();

	std::cout << "Created Illuminance Volume at: " << filePath << std::endl;
}
