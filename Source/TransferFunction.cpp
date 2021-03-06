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

#include "TransferFunction.h"

QTransferFunction gTransferFunction;

// Compare two transfer function nodes by intensity
bool CompareNodes(QNode NodeA, QNode NodeB)
{
	return NodeA.GetIntensity() < NodeB.GetIntensity();
}

QTransferFunction::QTransferFunction(QObject* pParent, const QString& Name) :
	QPresetXML(pParent),
	m_Nodes(),
	m_pSelectedNode(NULL),
	m_DensityScale(5.0f),
	m_ShadingType(1),
	m_GradientFactor(0.0f),
	m_AlgorithmType(0),
	m_ScatterType(2),
	m_NrOfBounces(1),
	m_PostProcessingSteps(31),
	m_PrimaryStepSize(3.0f),
	m_SecondaryStepSize(3.0f),
	m_ScatteringHeadstart(0),
	m_MakeFloodFill(false),
	m_OpacityWeight(100),
	m_DirectionWeight(1)
{
}

QTransferFunction::QTransferFunction(const QTransferFunction& Other)
{
	*this = Other;
};

QTransferFunction& QTransferFunction::operator = (const QTransferFunction& Other)			
{
	QPresetXML::operator=(Other);

	blockSignals(true);
	
	m_Nodes			= Other.m_Nodes;
	m_pSelectedNode	= Other.m_pSelectedNode;

	// Notify us when the nodes change
	for (int i = 0; i < m_Nodes.size(); i++)
		connect(&m_Nodes[i], SIGNAL(NodeChanged(QNode*)), this, SLOT(OnNodeChanged(QNode*)));

	m_DensityScale		= Other.m_DensityScale;
	m_ShadingType		= Other.m_ShadingType;
	m_GradientFactor	= Other.m_GradientFactor;
	m_AlgorithmType		= Other.m_AlgorithmType;
	m_ScatterType		= Other.m_ScatterType;
	m_NrOfBounces		= Other.m_NrOfBounces;
	m_PostProcessingSteps = Other.m_PostProcessingSteps;
	m_PrimaryStepSize	= Other.m_PrimaryStepSize;
	m_SecondaryStepSize = Other.m_SecondaryStepSize;
	m_ScatteringHeadstart = Other.m_ScatteringHeadstart;
	m_OpacityWeight = Other.m_OpacityWeight;
	m_DirectionWeight = Other.m_DirectionWeight;
	m_MakeFloodFill = Other.m_MakeFloodFill;

	// Update node's range
	UpdateNodeRanges();

	blockSignals(false);

	// Notify others that the function has changed selection has changed
	emit FunctionChanged();
	emit SettingsChanged();

	SetSelectedNode(NULL);

	return *this;
}

void QTransferFunction::OnNodeChanged(QNode* pNode)
{
	// Update node's range
	UpdateNodeRanges();

	emit FunctionChanged();

	SetDirty();
}

void QTransferFunction::SetSelectedNode(QNode* pSelectedNode)
{
	m_pSelectedNode = pSelectedNode;
	emit SelectionChanged(m_pSelectedNode);
}

void QTransferFunction::SetSelectedNode(const int& Index)
{
	if (m_Nodes.size() <= 0)
	{
		m_pSelectedNode = NULL;
	}
	else
	{
		// Compute new index
		const int NewIndex = qMin(m_Nodes.size() - 1, qMax(0, Index));

		// Set selected node
		m_pSelectedNode = &m_Nodes[NewIndex];
	}

	// Notify others that our selection has changed
	emit SelectionChanged(m_pSelectedNode);
}

QNode* QTransferFunction::GetSelectedNode(void)
{
	return m_pSelectedNode;
}

void QTransferFunction::SelectFirstNode(void)
{
	if (m_Nodes.size() == 0)
		return;

	SetSelectedNode(&m_Nodes[0]);
}

void QTransferFunction::SelectPreviousNode(void)
{
	if (!m_pSelectedNode)
		return;

	int Index = m_Nodes.indexOf(*GetSelectedNode());

	if (Index < 0)
		return;

	// Compute new index
	const int NewIndex = qMin(m_Nodes.size() - 1, qMax(0, Index - 1));

	// Set selected node
	SetSelectedNode(&m_Nodes[NewIndex]);
}

void QTransferFunction::SelectNextNode(void)
{
	if (!m_pSelectedNode)
		return;

	int Index = m_Nodes.indexOf(*GetSelectedNode());

	if (Index < 0)
		return;

	// Compute new index
	const int NewIndex = qMin(m_Nodes.size() - 1, qMax(0, Index + 1));

	// Set selected node
	SetSelectedNode(&m_Nodes[NewIndex]);
}

void QTransferFunction::SelectLastNode(void)
{
	if (m_Nodes.size() == 0)
		return;

	SetSelectedNode(&m_Nodes[m_Nodes.size() - 1]);
}

int	QTransferFunction::GetNodeIndex(QNode* pNode)
{
	if (pNode == NULL)
		return -1;
	
	return m_Nodes.indexOf(*pNode);
}

void QTransferFunction::AddNode(const float& Intensity, const float& Opacity, const QColor& Diffuse, const QColor& Specular, const QColor& Emission, const float& Roughness)
{
	AddNode(QNode(this, Intensity, Opacity, Diffuse, Specular, Emission, Roughness));
}

void QTransferFunction::AddNode(const QNode& Node)
{
	// Add the node to the list
	m_Nodes.append(Node);

	// Cache node
	QNode& CacheNode = m_Nodes.back();

	// Sort the transfer function nodes based on intensity
	qSort(m_Nodes.begin(), m_Nodes.end(), CompareNodes);

	// Update ID's
	for (int i = 0; i < m_Nodes.size(); i++)
		m_Nodes[i].m_ID = i;

	// Update ranges
	UpdateNodeRanges();

	// Notify us when the node changes
	connect(&CacheNode, SIGNAL(NodeChanged(QNode*)), this, SLOT(OnNodeChanged(QNode*)));

	for (int i = 0; i < m_Nodes.size(); i++)
	{
		if (Node.GetIntensity() == m_Nodes[i].GetIntensity())
			SetSelectedNode(&m_Nodes[i]);
	}

	// Inform others that the transfer function has changed
	emit FunctionChanged();

	if (!signalsBlocked())
		Log("Inserted node", "layer-select-point");
}

void QTransferFunction::RemoveNode(QNode* pNode)
{
	if (!pNode)
		return;

	// Remove the connection
	disconnect(pNode, SIGNAL(NodeChanged(QNode*)), this, SLOT(OnNodeChanged(QNode*)));

	// Node index of the to be removed node
	int NodeIndex = m_Nodes.indexOf(*pNode);

	// Remove from list and memory
	m_Nodes.removeOne(*pNode);

	// Update ID's
	for (int i = 0; i < m_Nodes.size(); i++)
		m_Nodes[i].m_ID = i;

	// Update node's range
	UpdateNodeRanges();

	// Select the previous node
	NodeIndex = qMax(0, NodeIndex - 1);

	SetSelectedNode(NodeIndex);

	// Inform others that the transfer function has changed
	emit FunctionChanged();

	Log("Removed node", "layer-select-point");
}

void QTransferFunction::UpdateNodeRanges(void)
{
	// Compute the node ranges
	for (int i = 0; i < m_Nodes.size(); i++)
	{
		QNode& Node = m_Nodes[i];

		if (i == 0)
		{
			Node.SetMinX(0.0f);
			Node.SetMaxX(0.0f);
		}
		else if (i == (m_Nodes.size() - 1))
		{
			Node.SetMinX(1.0f);
			Node.SetMaxX(1.0f);
		}
		else
		{
			QNode& NodeLeft		= m_Nodes[i - 1];
			QNode& NodeRight	= m_Nodes[i + 1];

			Node.SetMinX(NodeLeft.GetIntensity());
			Node.SetMaxX(NodeRight.GetIntensity());
		}
	}
}

const QNodeList& QTransferFunction::GetNodes(void) const
{
	return m_Nodes;
}

QNode& QTransferFunction::GetNode(const int& Index)
{
	return m_Nodes[Index];
}

float QTransferFunction::GetDensityScale(void) const
{
	return m_DensityScale;
}

void QTransferFunction::SetDensityScale(const float& DensityScale)
{
	if (DensityScale == m_DensityScale)
		return;

	m_DensityScale = DensityScale;

	emit FunctionChanged();
}

int QTransferFunction::GetShadingType(void) const
{
	return m_ShadingType;
}

void QTransferFunction::SetShadingType(const int& ShadingType)
{
	if (ShadingType == m_ShadingType)
		return;

	m_ShadingType = ShadingType;

	emit SettingsChanged();
}

float QTransferFunction::GetGradientFactor(void) const
{
	return m_GradientFactor;
}

void QTransferFunction::SetGradientFactor(const float& GradientFactor)
{
	if (GradientFactor == m_GradientFactor)
		return;

	m_GradientFactor = GradientFactor;

	emit FunctionChanged();
}

int QTransferFunction::GetAlgorithmType(void) const
{
	return m_AlgorithmType;
}

void QTransferFunction::SetAlgorithmType(const int& AlgorithmType)
{
	if (AlgorithmType == m_AlgorithmType)
		return;

	m_AlgorithmType = AlgorithmType;

	emit SettingsChanged();
}

int QTransferFunction::GetScatterType(void) const
{
	return m_ScatterType;
}

void QTransferFunction::SetScatterType(const int& ScatterType)
{
	if (ScatterType == m_ScatterType)
		return;

	m_ScatterType = ScatterType;

	emit SettingsChanged();
}

int QTransferFunction::GetNrOfBounces(void) const
{
	return m_NrOfBounces;
}

void QTransferFunction::SetNrOfBounces(const int& NrOfBounces)
{
	if (NrOfBounces == m_NrOfBounces)
		return;

	m_NrOfBounces = NrOfBounces;

	emit SettingsChanged();
}

short QTransferFunction::GetPostProcessingSteps(void) const
{
	return m_PostProcessingSteps;
}

void QTransferFunction::SetPostProcessingSteps(const short& PostProcessingSteps)
{
	if (PostProcessingSteps == m_PostProcessingSteps)
		return;

	m_PostProcessingSteps = PostProcessingSteps;

	emit SettingsChanged();
}

float QTransferFunction::GetPrimaryStepSize(void) const
{
	return m_PrimaryStepSize;
}

void QTransferFunction::SetPrimaryStepSize(const float& PrimaryStepSize)
{
	if (PrimaryStepSize == m_PrimaryStepSize)
		return;

	m_PrimaryStepSize = PrimaryStepSize;

	emit SettingsChanged();
}

float QTransferFunction::GetSecondarStepSize(void) const
{
	return m_SecondaryStepSize;
}

void QTransferFunction::SetSecondaryStepSize(const float& SecondaryStepSize)
{
	if (SecondaryStepSize == m_SecondaryStepSize)
		return;

	m_SecondaryStepSize = SecondaryStepSize;

	emit SettingsChanged();
}

float QTransferFunction::GetScatteringHeadstart(void) const
{
	return m_ScatteringHeadstart;
}

void QTransferFunction::SetScatteringHeadstart(const float& ScatteringHeadstart)
{
	if (ScatteringHeadstart == m_ScatteringHeadstart)
		return;

	m_ScatteringHeadstart = ScatteringHeadstart;

	emit SettingsChanged();
}

float QTransferFunction::GetOpacityWeight(void) const
{
	return m_OpacityWeight;
}

void QTransferFunction::SetOpacityWeight(const float& OpacityWeight)
{
	if (OpacityWeight == m_OpacityWeight)
		return;

	m_OpacityWeight = OpacityWeight;

	//emit SettingsChanged();
}

float QTransferFunction::GetDirectionWeight(void) const
{
	return m_DirectionWeight;
}

void QTransferFunction::SetDirectionWeight(const float& DirectionWeight)
{
	if (DirectionWeight == m_DirectionWeight)
		return;

	m_DirectionWeight = DirectionWeight;

	//emit SettingsChanged();
}

bool QTransferFunction::GetMakeFloodFill(void) const
{
	return m_MakeFloodFill;
}

void QTransferFunction::setMakeFloodFill(const bool& MakeFloodFill)
{
	m_MakeFloodFill = MakeFloodFill;

	if (m_MakeFloodFill)
		emit SettingsChanged();
}


void QTransferFunction::ReadXML(QDomElement& Parent)
{
	QPresetXML::ReadXML(Parent);

	QDomElement Nodes = Parent.firstChild().toElement();

	blockSignals(true);

	// Read child nodes
	for (QDomNode DomNode = Nodes.firstChild(); !DomNode.isNull(); DomNode = DomNode.nextSibling())
	{
		// Create new node
		QNode Node(this);

		// Load preset into it
		Node.ReadXML(DomNode.toElement());

		// Add the node to the list
		AddNode(Node);
	}

	UpdateNodeRanges();
	
	m_DensityScale		= Parent.firstChildElement("DensityScale").attribute("Value").toFloat();
	if (m_DensityScale == 0)
		m_DensityScale = 100;
	m_ShadingType		= Parent.firstChildElement("ShadingType").attribute("Value").toInt();
	m_GradientFactor	= Parent.firstChildElement("GradientFactor").attribute("Value").toFloat();
	m_AlgorithmType		= Parent.firstChildElement("AlgorithmType").attribute("Value").toInt();
	m_ScatterType		= Parent.firstChildElement("ScatterType").attribute("Value").toInt();
	m_NrOfBounces		= Parent.firstChildElement("NrOfBounces").attribute("Value").toInt();
	if (m_NrOfBounces == 0)
		m_NrOfBounces = 1;
	m_PostProcessingSteps = Parent.firstChildElement("PostProcessingSteps").attribute("Value").toShort();
	m_PrimaryStepSize = Parent.firstChildElement("PrimaryStepSize").attribute("Value").toFloat();
	if (m_PrimaryStepSize == 0)
		m_PrimaryStepSize = 3;
	m_SecondaryStepSize = Parent.firstChildElement("SecondaryStepSize").attribute("Value").toFloat();
	if (m_SecondaryStepSize == 0)
		m_SecondaryStepSize = 3;
	m_ScatteringHeadstart = Parent.firstChildElement("ScatteringHeadStart").attribute("Value").toFloat();

	blockSignals(false);

	// Inform others that the transfer function has changed
	emit FunctionChanged();
	emit SettingsChanged();
}

QDomElement QTransferFunction::WriteXML(QDomDocument& DOM, QDomElement& Parent)
{
	// Preset
	QDomElement Preset = DOM.createElement("Preset");
	Parent.appendChild(Preset);

	QPresetXML::WriteXML(DOM, Preset);

	Parent.appendChild(Preset);

	QDomElement Nodes = DOM.createElement("Nodes");
	Preset.appendChild(Nodes);

	for (int i = 0; i < m_Nodes.size(); i++)
		m_Nodes[i].WriteXML(DOM, Nodes);

	QDomElement DensityScale = DOM.createElement("DensityScale");
	DensityScale.setAttribute("Value", GetDensityScale());
	Preset.appendChild(DensityScale);

	QDomElement ShadingType = DOM.createElement("ShadingType");
	ShadingType.setAttribute("Value", GetShadingType());
	Preset.appendChild(ShadingType);

	QDomElement GradientFactor = DOM.createElement("GradientFactor");
	GradientFactor.setAttribute("Value", GetGradientFactor());
	Preset.appendChild(GradientFactor);

	QDomElement AlgorithmType = DOM.createElement("AlgorithmType");
	AlgorithmType.setAttribute("Value", GetAlgorithmType());
	Preset.appendChild(AlgorithmType);

	QDomElement ScatterType = DOM.createElement("ScatterType");
	ScatterType.setAttribute("Value", GetScatterType());
	Preset.appendChild(ScatterType);

	QDomElement NrOfBounces = DOM.createElement("NrOfBounces");
	NrOfBounces.setAttribute("Value", GetNrOfBounces());
	Preset.appendChild(NrOfBounces);
	
	QDomElement PostProcessingSteps = DOM.createElement("PostProcessingSteps");
	PostProcessingSteps.setAttribute("Value", GetPostProcessingSteps());
	Preset.appendChild(PostProcessingSteps);
	
	QDomElement PrimaryStepSize = DOM.createElement("PrimaryStepSize");
	PrimaryStepSize.setAttribute("Value", GetPrimaryStepSize());
	Preset.appendChild(PrimaryStepSize);

	QDomElement SecondaryStepSize = DOM.createElement("SecondaryStepSize");
	SecondaryStepSize.setAttribute("Value", GetSecondarStepSize());
	Preset.appendChild(SecondaryStepSize);
	
	QDomElement ScatteringHeadStart = DOM.createElement("ScatteringHeadStart");
	ScatteringHeadStart.setAttribute("Value", GetScatteringHeadstart());
	Preset.appendChild(ScatteringHeadStart);
	
	return Preset;
}

QTransferFunction QTransferFunction::Default(void)
{
	QTransferFunction DefaultTransferFunction;

	DefaultTransferFunction.SetName("Default");
	DefaultTransferFunction.AddNode(0.0f, 0.0f, Qt::gray, QColor(10, 10, 10), Qt::black, 1.0f);
	DefaultTransferFunction.AddNode(0.3f, 0.0f, Qt::gray, QColor(10, 10, 10), Qt::black, 1.0f);
	DefaultTransferFunction.AddNode(0.7f, 1.0f, Qt::gray, QColor(10, 10, 10), Qt::black, 1.0f);
	DefaultTransferFunction.AddNode(1.0f, 1.0f, Qt::gray, QColor(10, 10, 10), Qt::black, 1.0f);

	DefaultTransferFunction.SetDensityScale(100.0f);
	DefaultTransferFunction.SetShadingType(2);
	DefaultTransferFunction.SetGradientFactor(10.0f);
	DefaultTransferFunction.SetPostProcessingSteps(31);
	DefaultTransferFunction.SetPrimaryStepSize(3);
	DefaultTransferFunction.SetSecondaryStepSize(3);
	DefaultTransferFunction.SetNrOfBounces(1);
	
	return DefaultTransferFunction;
}