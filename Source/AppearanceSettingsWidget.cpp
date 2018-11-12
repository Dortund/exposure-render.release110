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

#include "AppearanceSettingsWidget.h"
#include "TransferFunction.h"
#include "RenderThread.h"
#include "Scene.h"

QAppearanceSettingsWidget::QAppearanceSettingsWidget(QWidget* pParent) :
	QGroupBox(pParent),
	m_MainLayout(),
	m_DensityScaleSlider(),
	m_DensityScaleSpinner(),
	m_ShadingType(),
	m_GradientFactorLabel(),
	m_GradientFactorSlider(),
	m_GradientFactorSpinner(),
	m_StepSizePrimaryRaySlider(),
	m_StepSizePrimaryRaySpinner(),
	m_StepSizeSecondaryRaySlider(),
	m_StepSizeSecondaryRaySpinner(),
	m_AlgorithmType(),
	m_DoBlur(),
	m_DoEstimate(),
	m_DoToneMap(),
	m_DoDenoise(),
	m_DoOffset()
{
	setLayout(&m_MainLayout);

	m_MainLayout.addWidget(new QLabel("Algorithm Type"), 1, 0);

	m_AlgorithmType.addItem("Single Scattering", 0);
	m_AlgorithmType.addItem("Multiple Scattering", 1);
	m_AlgorithmType.addItem("Pre-Random", 2);
	m_AlgorithmType.addItem("Opacity Gradient", 3);
	m_MainLayout.addWidget(&m_AlgorithmType, 1, 1, 1, 2);

	QObject::connect(&m_AlgorithmType, SIGNAL(currentIndexChanged(int)), this, SLOT(onSetAlgorithmType(int)));

	m_MainLayout.addWidget(new QLabel("Density Scale"), 2, 0);

	m_DensityScaleSlider.setOrientation(Qt::Horizontal);
	m_DensityScaleSlider.setRange(0.001, 100.0);
	m_DensityScaleSlider.setValue(1.0);
	m_MainLayout.addWidget(&m_DensityScaleSlider, 2, 1);

	m_DensityScaleSpinner.setRange(0.001, 100.0);
	m_DensityScaleSpinner.setDecimals(3);
	m_MainLayout.addWidget(&m_DensityScaleSpinner, 2, 2);

	m_MainLayout.addWidget(new QLabel("Shading Type"), 3, 0);

	m_ShadingType.addItem("BRDF Only", 0);
	m_ShadingType.addItem("Phase Function Only", 1);
	m_ShadingType.addItem("Hybrid", 2);
	m_MainLayout.addWidget(&m_ShadingType, 3, 1, 1, 2);

	m_GradientFactorLabel.setText("Gradient Factor");
	m_MainLayout.addWidget(&m_GradientFactorLabel, 4, 0);
	
	m_GradientFactorSlider.setRange(0.001, 100.0);
	m_GradientFactorSlider.setValue(100.0);

	m_MainLayout.addWidget(&m_GradientFactorSlider, 4, 1);

	m_GradientFactorSpinner.setRange(0.001, 100.0);
	m_GradientFactorSpinner.setDecimals(3);

	m_MainLayout.addWidget(&m_GradientFactorSpinner, 4, 2);

	QObject::connect(&m_DensityScaleSlider, SIGNAL(valueChanged(double)), &m_DensityScaleSpinner, SLOT(setValue(double)));
	QObject::connect(&m_DensityScaleSpinner, SIGNAL(valueChanged(double)), &m_DensityScaleSlider, SLOT(setValue(double)));
	QObject::connect(&m_DensityScaleSlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetDensityScale(double)));

	QObject::connect(&m_GradientFactorSlider, SIGNAL(valueChanged(double)), &m_GradientFactorSpinner, SLOT(setValue(double)));
	QObject::connect(&m_GradientFactorSpinner, SIGNAL(valueChanged(double)), &m_GradientFactorSlider, SLOT(setValue(double)));
	QObject::connect(&m_GradientFactorSlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetGradientFactor(double)));

	m_MainLayout.addWidget(new QLabel("Primary Step Size"), 5, 0);

	m_StepSizePrimaryRaySlider.setRange(1.0, 10.0);

	m_MainLayout.addWidget(&m_StepSizePrimaryRaySlider, 5, 1);

	m_StepSizePrimaryRaySpinner.setRange(1.0, 10.0);
	m_StepSizePrimaryRaySpinner.setDecimals(2);

	m_MainLayout.addWidget(&m_StepSizePrimaryRaySpinner, 5, 2);

	QObject::connect(&m_StepSizePrimaryRaySlider, SIGNAL(valueChanged(double)), &m_StepSizePrimaryRaySpinner, SLOT(setValue(double)));
	QObject::connect(&m_StepSizePrimaryRaySpinner, SIGNAL(valueChanged(double)), &m_StepSizePrimaryRaySlider, SLOT(setValue(double)));
	QObject::connect(&m_StepSizePrimaryRaySlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetStepSizePrimaryRay(double)));

	m_MainLayout.addWidget(new QLabel("Secondary Step Size"), 6, 0);

	m_StepSizeSecondaryRaySlider.setRange(1.0, 10.0);

	m_MainLayout.addWidget(&m_StepSizeSecondaryRaySlider, 6, 1);

	m_StepSizeSecondaryRaySpinner.setRange(1.0, 10.0);
	m_StepSizeSecondaryRaySpinner.setDecimals(2);

	m_MainLayout.addWidget(&m_StepSizeSecondaryRaySpinner, 6, 2);

	QObject::connect(&m_StepSizeSecondaryRaySlider, SIGNAL(valueChanged(double)), &m_StepSizeSecondaryRaySpinner, SLOT(setValue(double)));
	QObject::connect(&m_StepSizeSecondaryRaySpinner, SIGNAL(valueChanged(double)), &m_StepSizeSecondaryRaySlider, SLOT(setValue(double)));
	QObject::connect(&m_StepSizeSecondaryRaySlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetStepSizeSecondaryRay(double)));


	QObject::connect(&m_ShadingType, SIGNAL(currentIndexChanged(int)), this, SLOT(OnSetShadingType(int)));
	QObject::connect(&gStatus, SIGNAL(RenderBegin()), this, SLOT(OnRenderBegin()));
	QObject::connect(&gTransferFunction, SIGNAL(Changed()), this, SLOT(OnTransferFunctionChanged()));

	m_MainLayout.addWidget(new QLabel("Random Offset"), 7, 0);

	m_DoOffset.setChecked(true);
	m_MainLayout.addWidget(&m_DoOffset, 7, 1);

	QObject::connect(&m_DoOffset, SIGNAL(stateChanged(int)), this, SLOT(onDoOffsetChanged(int)));

	m_MainLayout.addWidget(new QLabel("Blur"), 8, 0);

	m_DoBlur.setChecked(true);
	m_MainLayout.addWidget(&m_DoBlur, 8, 1);

	QObject::connect(&m_DoBlur, SIGNAL(stateChanged(int)), this, SLOT(OnDoBlurChanged(int)));


	m_MainLayout.addWidget(new QLabel("Esitmate"), 9, 0);

	m_DoEstimate.setChecked(true);
	m_MainLayout.addWidget(&m_DoEstimate, 9, 1);

	QObject::connect(&m_DoEstimate, SIGNAL(stateChanged(int)), this, SLOT(OnDoEstimateChanged(int)));

	// TODO atm we always do tonemap. remove related code or make optional
	//m_MainLayout.addWidget(new QLabel("ToneMap"), 9, 0);

	m_DoToneMap.setChecked(true);
	//m_MainLayout.addWidget(&m_DoToneMap, 9, 1);

	QObject::connect(&m_DoToneMap, SIGNAL(stateChanged(int)), this, SLOT(OnDoToneMapChanged(int)));


	m_MainLayout.addWidget(new QLabel("Denoise"), 10, 0);

	m_DoDenoise.setChecked(true);
	m_MainLayout.addWidget(&m_DoDenoise, 10, 1);

	QObject::connect(&m_DoDenoise, SIGNAL(stateChanged(int)), this, SLOT(OnDoDenoiseChanged(int)));
}

void QAppearanceSettingsWidget::OnRenderBegin(void)
{
	m_DensityScaleSlider.setValue(gTransferFunction.GetDensityScale());
	m_ShadingType.setCurrentIndex(gTransferFunction.GetShadingType());
	m_GradientFactorSlider.setValue(gScene.m_GradientFactor);

	m_StepSizePrimaryRaySlider.setValue(gScene.m_StepSizeFactor, true);
	m_StepSizePrimaryRaySpinner.setValue(gScene.m_StepSizeFactor, true);
	m_StepSizeSecondaryRaySlider.setValue(gScene.m_StepSizeFactorShadow, true);
	m_StepSizeSecondaryRaySpinner.setValue(gScene.m_StepSizeFactorShadow, true);
}

void QAppearanceSettingsWidget::OnSetDensityScale(double DensityScale)
{
	gTransferFunction.SetDensityScale(DensityScale);
}

void QAppearanceSettingsWidget::OnSetShadingType(int Index)
{
	gTransferFunction.SetShadingType(Index);
	m_GradientFactorLabel.setEnabled(Index == 2);
	m_GradientFactorSlider.setEnabled(Index == 2);
	m_GradientFactorSpinner.setEnabled(Index == 2);
}

void QAppearanceSettingsWidget::onSetAlgorithmType(int Index)
{
	gTransferFunction.SetAlgorithmType(Index);
}

void QAppearanceSettingsWidget::OnSetGradientFactor(double GradientFactor)
{
	gTransferFunction.SetGradientFactor(GradientFactor);
}

void QAppearanceSettingsWidget::OnSetStepSizePrimaryRay(const double& StepSizePrimaryRay)
{
	gScene.m_StepSizeFactor = (float)StepSizePrimaryRay;
	gScene.m_DirtyFlags.SetFlag(RenderParamsDirty);
}

void QAppearanceSettingsWidget::OnSetStepSizeSecondaryRay(const double& StepSizeSecondaryRay)
{
	gScene.m_StepSizeFactorShadow = (float)StepSizeSecondaryRay;
	gScene.m_DirtyFlags.SetFlag(RenderParamsDirty);
}

void QAppearanceSettingsWidget::OnTransferFunctionChanged(void)
{
	m_DensityScaleSlider.setValue(gTransferFunction.GetDensityScale(), true);
	m_DensityScaleSlider.setValue(gTransferFunction.GetDensityScale(), true);
	m_ShadingType.setCurrentIndex(gTransferFunction.GetShadingType());
	m_GradientFactorSlider.setValue(gTransferFunction.GetGradientFactor(), true);
	m_GradientFactorSpinner.setValue(gTransferFunction.GetGradientFactor(), true);
}

void QAppearanceSettingsWidget::OnDoBlurChanged(int doBlur) {
	gScene.m_PostProcessingSteps ^= 1;
	gScene.m_DirtyFlags.SetFlag(RenderParamsDirty);
}

void QAppearanceSettingsWidget::OnDoEstimateChanged(int doEstimate) {
	gScene.m_PostProcessingSteps ^= 2;
	gScene.m_DirtyFlags.SetFlag(RenderParamsDirty);
}

void QAppearanceSettingsWidget::OnDoToneMapChanged(int doToneMap) {
	gScene.m_PostProcessingSteps ^= 4;
	gScene.m_DirtyFlags.SetFlag(RenderParamsDirty);
}

void QAppearanceSettingsWidget::OnDoDenoiseChanged(int doDenoise) {
	gScene.m_PostProcessingSteps ^= 8;
	gScene.m_DirtyFlags.SetFlag(RenderParamsDirty);
}

void QAppearanceSettingsWidget::onDoOffsetChanged(int doOffset) {
	gScene.m_PostProcessingSteps ^= 16;
	gScene.m_DirtyFlags.SetFlag(RenderParamsDirty);
}