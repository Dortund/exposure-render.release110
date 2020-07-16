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
	m_DoOffset(),
	m_MaxBouncesSpinner(),
	m_MakeFloodFillButton("Make New Floodfill")
{
	setLayout(&m_MainLayout);

	int i = 0;

	m_MainLayout.addWidget(new QLabel("Algorithm Type"), i, 0);

	m_AlgorithmType.addItem("Single Scattering", 0);
	m_AlgorithmType.addItem("Multiple Scattering", 1);
	m_AlgorithmType.addItem("Pre-Random", 2);
	m_AlgorithmType.addItem("Opacity Gradient", 3);
	m_AlgorithmType.addItem("Flood Fill", 4);
	m_AlgorithmType.addItem("Nr of Bounces", 5);
	m_AlgorithmType.addItem("Final Throughput", 6);
	m_AlgorithmType.addItem("Bounce Type Order", 7);
	m_AlgorithmType.addItem("Property Based", 8);
	m_AlgorithmType.addItem("Property Based Normal", 9);
	m_AlgorithmType.addItem("Property Based Opacity/roughness/magnitude", 10);
	m_AlgorithmType.addItem("Property Based Fractions", 11);
	m_AlgorithmType.addItem("Property Based Voxels", 12);
	m_AlgorithmType.addItem("Property Based Intensity", 13);
	m_AlgorithmType.addItem("Property Based Coordinate", 14);
	m_AlgorithmType.addItem("Multi Coordinate", 15);
	m_AlgorithmType.addItem("Travel Distance", 16);
	m_AlgorithmType.addItem("Likeliest Light Direction", 17);
	m_AlgorithmType.addItem("Likeliest Light Direction Cube", 18);
	m_AlgorithmType.addItem("Light Path Values", 19);
	m_AlgorithmType.addItem("Chosen Light Path", 20);
	m_AlgorithmType.addItem("Light Path Variance", 21);
	m_AlgorithmType.addItem("Property Based Final Throughput", 22);
	m_MainLayout.addWidget(&m_AlgorithmType, i++, 1, 1, 2);

	QObject::connect(&m_AlgorithmType, SIGNAL(currentIndexChanged(int)), this, SLOT(OnSetAlgorithmType(int)));

	m_MainLayout.addWidget(new QLabel("Shading Type"), i, 0);

	m_ShadingType.addItem("BRDF Only", 0);
	m_ShadingType.addItem("Phase Function Only", 1);
	m_ShadingType.addItem("Hybrid", 2);
	m_ShadingType.addItem("Light Paths", 3);
	m_ShadingType.addItem("Light Paths Octo", 4);
	m_ShadingType.addItem("Light Paths Octo Gradient", 5);
	m_MainLayout.addWidget(&m_ShadingType, i++, 1, 1, 2);

	QObject::connect(&m_ShadingType, SIGNAL(currentIndexChanged(int)), this, SLOT(OnSetShadingType(int)));

	m_MainLayout.addWidget(new QLabel("Scattering Type"), i, 0);

	m_ScatterType.addItem("BRDF Only", 0);
	m_ScatterType.addItem("Phase Function Only", 1);
	m_ScatterType.addItem("Hybrid", 2);
	m_ScatterType.addItem("Light Paths", 3);
	m_ScatterType.addItem("Light Paths Octo", 4);
	m_ScatterType.addItem("Light Paths Octo Gradient", 5);
	m_ScatterType.setCurrentIndex(2);
	m_MainLayout.addWidget(&m_ScatterType, i++, 1, 1, 2);

	QObject::connect(&m_ScatterType, SIGNAL(currentIndexChanged(int)), this, SLOT(OnSetScatteringType(int)));

	QLabel* lblDensityScale = new QLabel("Density Scale");
	lblDensityScale->setToolTip("Higher values virtually scale up the opacity and reduce the overal distance a ray can travel through the volume before colliding");

	m_MainLayout.addWidget(lblDensityScale, i, 0);

	m_DensityScaleSlider.setOrientation(Qt::Horizontal);
	m_DensityScaleSlider.setRange(0.001, 1000.0);
	m_DensityScaleSlider.setValue(1.0);
	m_MainLayout.addWidget(&m_DensityScaleSlider, i, 1);

	m_DensityScaleSpinner.setRange(0.001, 1000.0);
	m_DensityScaleSpinner.setDecimals(3);
	m_MainLayout.addWidget(&m_DensityScaleSpinner, i++, 2);

	QObject::connect(&m_DensityScaleSlider, SIGNAL(valueChanged(double)), &m_DensityScaleSpinner, SLOT(setValue(double)));
	QObject::connect(&m_DensityScaleSpinner, SIGNAL(valueChanged(double)), &m_DensityScaleSlider, SLOT(setValue(double)));
	QObject::connect(&m_DensityScaleSlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetDensityScale(double)));

	m_GradientFactorLabel.setText("Gradient Factor");
	m_GradientFactorLabel.setToolTip("Multiplies the gradient with this factor when calculating getting the gradient for hybrid shading");
	m_MainLayout.addWidget(&m_GradientFactorLabel, i, 0);
	
	m_GradientFactorSlider.setRange(0.001, 100.0);
	m_GradientFactorSlider.setValue(100.0);

	m_MainLayout.addWidget(&m_GradientFactorSlider, i, 1);

	m_GradientFactorSpinner.setRange(0.001, 100.0);
	m_GradientFactorSpinner.setDecimals(3);

	m_MainLayout.addWidget(&m_GradientFactorSpinner, i++, 2);

	QObject::connect(&m_GradientFactorSlider, SIGNAL(valueChanged(double)), &m_GradientFactorSpinner, SLOT(setValue(double)));
	QObject::connect(&m_GradientFactorSpinner, SIGNAL(valueChanged(double)), &m_GradientFactorSlider, SLOT(setValue(double)));
	QObject::connect(&m_GradientFactorSlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetGradientFactor(double)));

	QLabel* lblPrimStepSize = new QLabel("Primary Step Size");
	lblPrimStepSize->setToolTip("Sets the stepsize in voxels for primary rays");
	m_MainLayout.addWidget(lblPrimStepSize, i, 0);

	const double stepMin = 0.01;

	m_StepSizePrimaryRaySlider.setRange(stepMin, 10.0);

	m_MainLayout.addWidget(&m_StepSizePrimaryRaySlider, i, 1);

	m_StepSizePrimaryRaySpinner.setRange(stepMin, 10.0);
	m_StepSizePrimaryRaySpinner.setDecimals(2);

	m_MainLayout.addWidget(&m_StepSizePrimaryRaySpinner, i++, 2);

	QObject::connect(&m_StepSizePrimaryRaySlider, SIGNAL(valueChanged(double)), &m_StepSizePrimaryRaySpinner, SLOT(setValue(double)));
	QObject::connect(&m_StepSizePrimaryRaySpinner, SIGNAL(valueChanged(double)), &m_StepSizePrimaryRaySlider, SLOT(setValue(double)));
	QObject::connect(&m_StepSizePrimaryRaySlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetStepSizePrimaryRay(double)));

	QLabel* lblSecStepSize = new QLabel("Secondary Step Size");
	lblSecStepSize->setToolTip("Sets the stepsize in voxels for secondary rays");
	m_MainLayout.addWidget(lblSecStepSize, i, 0);

	m_StepSizeSecondaryRaySlider.setRange(stepMin, 10.0);

	m_MainLayout.addWidget(&m_StepSizeSecondaryRaySlider, i, 1);

	m_StepSizeSecondaryRaySpinner.setRange(stepMin, 10.0);
	m_StepSizeSecondaryRaySpinner.setDecimals(2);

	m_MainLayout.addWidget(&m_StepSizeSecondaryRaySpinner, i++, 2);

	QObject::connect(&m_StepSizeSecondaryRaySlider, SIGNAL(valueChanged(double)), &m_StepSizeSecondaryRaySpinner, SLOT(setValue(double)));
	QObject::connect(&m_StepSizeSecondaryRaySpinner, SIGNAL(valueChanged(double)), &m_StepSizeSecondaryRaySlider, SLOT(setValue(double)));
	QObject::connect(&m_StepSizeSecondaryRaySlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetStepSizeSecondaryRay(double)));

	QLabel* lblScatHeadstart = new QLabel("Scattering Ray Headstart");
	lblScatHeadstart->setToolTip("Sets the raymarching starting distance from a scattering point");
	m_MainLayout.addWidget(lblScatHeadstart, i, 0);

	m_ScatteringHeadstartSlider.setRange(0, 10.0);

	m_MainLayout.addWidget(&m_ScatteringHeadstartSlider, i, 1);

	m_ScatteringHeadstartSpinner.setRange(0, 10.0);
	m_ScatteringHeadstartSpinner.setDecimals(2);

	m_MainLayout.addWidget(&m_ScatteringHeadstartSpinner, i++, 2);

	QObject::connect(&m_ScatteringHeadstartSlider, SIGNAL(valueChanged(double)), &m_ScatteringHeadstartSpinner, SLOT(setValue(double)));
	QObject::connect(&m_ScatteringHeadstartSpinner, SIGNAL(valueChanged(double)), &m_ScatteringHeadstartSlider, SLOT(setValue(double)));
	QObject::connect(&m_ScatteringHeadstartSlider, SIGNAL(valueChanged(double)), this, SLOT(OnSetScatteringHeadstart(double)));

	m_MainLayout.addWidget(new QLabel("Random Offset"), i, 0);

	m_DoOffset.setChecked(true);
	m_MainLayout.addWidget(&m_DoOffset, i++, 1);

	QObject::connect(&m_DoOffset, SIGNAL(stateChanged(int)), this, SLOT(OnDoOffsetChanged(int)));

	m_MainLayout.addWidget(new QLabel("Blur"), i, 0);

	m_DoBlur.setChecked(true);
	m_MainLayout.addWidget(&m_DoBlur, i++, 1);

	QObject::connect(&m_DoBlur, SIGNAL(stateChanged(int)), this, SLOT(OnDoBlurChanged(int)));


	m_MainLayout.addWidget(new QLabel("Estimate"), i, 0);

	m_DoEstimate.setChecked(true);
	m_MainLayout.addWidget(&m_DoEstimate, i++, 1);

	QObject::connect(&m_DoEstimate, SIGNAL(stateChanged(int)), this, SLOT(OnDoEstimateChanged(int)));

	
	m_MainLayout.addWidget(new QLabel("ToneMap"), i, 0);

	m_DoToneMap.setChecked(true);
	m_MainLayout.addWidget(&m_DoToneMap, i++, 1);

	QObject::connect(&m_DoToneMap, SIGNAL(stateChanged(int)), this, SLOT(OnDoToneMapChanged(int)));


	m_MainLayout.addWidget(new QLabel("Denoise"), i, 0);

	m_DoDenoise.setChecked(true);
	m_MainLayout.addWidget(&m_DoDenoise, i++, 1);

	QObject::connect(&m_DoDenoise, SIGNAL(stateChanged(int)), this, SLOT(OnDoDenoiseChanged(int)));


	m_MainLayout.addWidget(new QLabel("max # of Bounces"), i, 0);

	m_MaxBouncesSpinner.setRange(1, 100);
	m_MaxBouncesSpinner.setDecimals(0);

	m_MainLayout.addWidget(&m_MaxBouncesSpinner, i++, 1);

	QObject::connect(&m_MaxBouncesSpinner, SIGNAL(valueChanged(double)), this, SLOT(OnSetMaxBounces(double)));

	QLabel* lblOpacityWeight = new QLabel("FF Opacity Weight");
	//lblOpacityWeight->setToolTip("Sets the raymarching starting distance from a scattering point");
	m_MainLayout.addWidget(lblOpacityWeight, i, 0);

	//m_ScatteringHeadstartSlider.setRange(0, 10.0);

	//m_MainLayout.addWidget(&m_ScatteringHeadstartSlider, i, 1);

	m_OpacityWeightSpinner.setRange(0, 100000000.0);
	m_OpacityWeightSpinner.setDecimals(0);

	m_MainLayout.addWidget(&m_OpacityWeightSpinner, i++, 2);

	//QObject::connect(&m_ScatteringHeadstartSlider, SIGNAL(valueChanged(double)), &m_ScatteringHeadstartSpinner, SLOT(setValue(double)));
	//QObject::connect(&m_ScatteringHeadstartSpinner, SIGNAL(valueChanged(double)), &m_ScatteringHeadstartSlider, SLOT(setValue(double)));
	QObject::connect(&m_OpacityWeightSpinner, SIGNAL(valueChanged(double)), this, SLOT(OnSetOpacityWeight(double)));

	m_MainLayout.addWidget(&m_MakeFloodFillButton, i, 0);
	QObject::connect(&m_MakeFloodFillButton, SIGNAL(released()), this, SLOT(OnMakeFloodfill()));


	QObject::connect(&gStatus, SIGNAL(RenderBegin()), this, SLOT(OnRenderBegin()));
	QObject::connect(&gTransferFunction, SIGNAL(FunctionChanged()), this, SLOT(OnTransferFunctionChanged()));
	QObject::connect(&gTransferFunction, SIGNAL(SettingsChanged()), this, SLOT(OnTransferFunctionSettingsChanged()));

}

void QAppearanceSettingsWidget::OnRenderBegin(void)
{
	OnTransferFunctionSettingsChanged();
	OnTransferFunctionChanged();
}

void QAppearanceSettingsWidget::OnTransferFunctionChanged(void)
{
	m_DensityScaleSlider.setValue(gTransferFunction.GetDensityScale(), true);
	m_DensityScaleSpinner.setValue(gTransferFunction.GetDensityScale(), true);
	m_GradientFactorSlider.setValue(gTransferFunction.GetGradientFactor(), true);
	m_GradientFactorSpinner.setValue(gTransferFunction.GetGradientFactor(), true);
}

void QAppearanceSettingsWidget::OnTransferFunctionSettingsChanged(void) {
	if (m_ShadingType.currentIndex() != gTransferFunction.GetShadingType())
		m_ShadingType.setCurrentIndex(gTransferFunction.GetShadingType());
	if (m_AlgorithmType.currentIndex() != gTransferFunction.GetAlgorithmType())
		m_AlgorithmType.setCurrentIndex(gTransferFunction.GetAlgorithmType());
	if (m_ScatterType.currentIndex() != gTransferFunction.GetScatterType())
		m_ScatterType.setCurrentIndex(gTransferFunction.GetScatterType());
	if (m_MaxBouncesSpinner.value() != gTransferFunction.GetNrOfBounces(), true)
		m_MaxBouncesSpinner.setValue(gTransferFunction.GetNrOfBounces(), true);
	if (m_DoBlur.isChecked() != (bool)(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::BLUR)) {
		m_DoBlur.blockSignals(true);
		m_DoBlur.setChecked(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::BLUR);
		m_DoBlur.blockSignals(false);
	}
	if (m_DoEstimate.isChecked() != (bool)(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::ESTIMATE)) {
		m_DoEstimate.blockSignals(true);
		m_DoEstimate.setChecked(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::ESTIMATE);
		m_DoEstimate.blockSignals(false);
	}
	if (m_DoToneMap.isChecked() != (bool)(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::TONE_MAP)) {
		m_DoToneMap.blockSignals(true);
		m_DoToneMap.setChecked(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::TONE_MAP);
		m_DoToneMap.blockSignals(false);
	}
	if (m_DoDenoise.isChecked() != (bool)(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::DENOISE)) {
		m_DoDenoise.blockSignals(true);
		m_DoDenoise.setChecked(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::DENOISE);
		m_DoDenoise.blockSignals(false);
	}
	if (m_DoOffset.isChecked() != (bool)(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::OFFSET)) {
		m_DoOffset.blockSignals(true);
		m_DoOffset.setChecked(gTransferFunction.GetPostProcessingSteps() & PostProcessingStepsEnum::OFFSET);
		m_DoOffset.blockSignals(false);
	}
	m_StepSizePrimaryRaySlider.setValue(gTransferFunction.GetPrimaryStepSize(), true);
	m_StepSizePrimaryRaySpinner.setValue(gTransferFunction.GetPrimaryStepSize(), true);
	m_StepSizeSecondaryRaySlider.setValue(gTransferFunction.GetSecondarStepSize(), true);
	m_StepSizeSecondaryRaySpinner.setValue(gTransferFunction.GetSecondarStepSize(), true);
	m_ScatteringHeadstartSlider.setValue(gTransferFunction.GetScatteringHeadstart(), true);
	m_ScatteringHeadstartSpinner.setValue(gTransferFunction.GetScatteringHeadstart(), true);
}

void QAppearanceSettingsWidget::OnMakeFloodfill() {
	gTransferFunction.setMakeFloodFill(true);
	m_MakeFloodFillButton.setEnabled(false);
}

void QAppearanceSettingsWidget::OnSetOpacityWeight(double OpacityWeight)
{
	gTransferFunction.SetOpacityWeight(OpacityWeight);
	m_MakeFloodFillButton.setEnabled(true);
}

void QAppearanceSettingsWidget::OnSetDirectionWeight(double DirectionWeight)
{
	gTransferFunction.SetDirectionWeight(DirectionWeight);
	m_MakeFloodFillButton.setEnabled(true);
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

void QAppearanceSettingsWidget::OnSetAlgorithmType(int Index)
{
	int val = m_AlgorithmType.itemData(Index).toInt();
	gTransferFunction.SetAlgorithmType(val);
}

void QAppearanceSettingsWidget::OnSetGradientFactor(double GradientFactor)
{
	gTransferFunction.SetGradientFactor(GradientFactor);
}

void QAppearanceSettingsWidget::OnSetStepSizePrimaryRay(const double& StepSizePrimaryRay)
{
	gTransferFunction.SetPrimaryStepSize((float)StepSizePrimaryRay);
}

void QAppearanceSettingsWidget::OnSetStepSizeSecondaryRay(const double& StepSizeSecondaryRay)
{
	gTransferFunction.SetSecondaryStepSize((float)StepSizeSecondaryRay);
}

void QAppearanceSettingsWidget::OnDoBlurChanged(int doBlur) {
	gTransferFunction.SetPostProcessingSteps(gTransferFunction.GetPostProcessingSteps() ^ PostProcessingStepsEnum::BLUR);
}

void QAppearanceSettingsWidget::OnDoEstimateChanged(int doEstimate) {
	gTransferFunction.SetPostProcessingSteps(gTransferFunction.GetPostProcessingSteps() ^ PostProcessingStepsEnum::ESTIMATE);
}

void QAppearanceSettingsWidget::OnDoToneMapChanged(int doToneMap) {
	gTransferFunction.SetPostProcessingSteps(gTransferFunction.GetPostProcessingSteps() ^ PostProcessingStepsEnum::TONE_MAP);
}

void QAppearanceSettingsWidget::OnDoDenoiseChanged(int doDenoise) {
	gTransferFunction.SetPostProcessingSteps(gTransferFunction.GetPostProcessingSteps() ^ PostProcessingStepsEnum::DENOISE);
}

void QAppearanceSettingsWidget::OnDoOffsetChanged(int doOffset) {
	gTransferFunction.SetPostProcessingSteps(gTransferFunction.GetPostProcessingSteps() ^ PostProcessingStepsEnum::OFFSET);
}

void QAppearanceSettingsWidget::OnSetMaxBounces(double nrOfBounces) {
	gTransferFunction.SetNrOfBounces((int)nrOfBounces);
}

void QAppearanceSettingsWidget::OnSetScatteringType(int index)
{
	gTransferFunction.SetScatterType(index);
}

void QAppearanceSettingsWidget::OnSetScatteringHeadstart(double ScatteringHeadstart) {
	gTransferFunction.SetScatteringHeadstart((float)ScatteringHeadstart);
}