/****************************************************************************
**
** Copyright (C) 2013 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the demonstration applications of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Digia.  For licensing terms and
** conditions see http://qt.digia.com/licensing.  For further information
** use the contact form at http://qt.digia.com/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** In addition, as a special exception, Digia gives you certain additional
** rights.  These rights are described in the Digia Qt LGPL Exception
** version 1.1, included in the file LGPL_EXCEPTION.txt in this package.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements will be
** met: http://www.gnu.org/copyleft/gpl.html.
**
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "controlwidget.h"

#include <QLayout>
//#include <QPainter>
//#include <QPainterPath>

const int alpha = 155;

ControlView::ControlView(QWidget *parent)
: QWidget(parent)
{
    setAttribute(Qt::WA_MouseTracking);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void ControlView::mousePressEvent(QMouseEvent *)
{
    //setDescriptionEnabled(false);
}

void ControlView::resizeEvent(QResizeEvent *e)
{
    //pts->setBoundingRect(rect());
    QWidget::resizeEvent(e);
}

#ifndef QT_NO_WHEELEVENT
void ControlView::wheelEvent(QWheelEvent *e)
{
    //m_scale += e->delta() / qreal(600);
    //m_scale = qMax(qreal(0.1), qMin(qreal(4), m_scale));
    //emit scaleChanged(int(m_scale*1000));
}
#endif

ControlWidget::ControlWidget(QWidget *parent)
    : QWidget(parent), textEditor(new QLineEdit)
{
    //setWindowTitle(tr("Affine Transformations"));

    /*view = new ControlView(this);
    view->setMinimumSize(0, 0);
    view->setMaximumSize(50, 100);*/

    this->setFixedWidth(200);

    QGroupBox *mainGroup = new QGroupBox(this);
    mainGroup->setFixedWidth(180);
    mainGroup->setTitle(tr("Control"));

    QGroupBox *rotateGroup = new QGroupBox(mainGroup);
    rotateGroup->setTitle(tr("Rotate"));
    QSlider *rotateSlider = new QSlider(Qt::Horizontal, rotateGroup);
    rotateSlider->setRange(0, 3600);
    rotateSlider->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

    QGroupBox *scaleGroup = new QGroupBox(mainGroup);
    scaleGroup->setTitle(tr("Scale"));
    QSlider *scaleSlider = new QSlider(Qt::Horizontal, scaleGroup);
    scaleSlider->setRange(1, 4000);
    scaleSlider->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

    QGroupBox *shearGroup = new QGroupBox(mainGroup);
    shearGroup->setTitle(tr("Shear"));
    QSlider *shearSlider = new QSlider(Qt::Horizontal, shearGroup);
    shearSlider->setRange(-990, 990);
    shearSlider->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

    QGroupBox *typeGroup = new QGroupBox(mainGroup);
    typeGroup->setTitle(tr("Type"));
    QRadioButton *vectorType = new QRadioButton(typeGroup);
    QRadioButton *pixmapType = new QRadioButton(typeGroup);
    QRadioButton *textType= new QRadioButton(typeGroup);
    vectorType->setText(tr("Vector Image"));
    pixmapType->setText(tr("Pixmap"));
    textType->setText(tr("Text"));

    QPushButton *resetButton = new QPushButton(mainGroup);
    resetButton->setText(tr("Reset Transform"));

    QPushButton *animateButton = new QPushButton(mainGroup);
    animateButton->setText(tr("Animate"));
    animateButton->setCheckable(true);

    QPushButton *showSourceButton = new QPushButton(mainGroup);
    showSourceButton->setText(tr("Show Source"));
#ifdef QT_OPENGL_SUPPORT
    QPushButton *enableOpenGLButton = new QPushButton(mainGroup);
    enableOpenGLButton->setText(tr("Use OpenGL"));
    enableOpenGLButton->setCheckable(true);
    //enableOpenGLButton->setChecked(view->usesOpenGL());
    //if (!QGLFormat::hasOpenGL())
    //    enableOpenGLButton->hide();
#endif
    QPushButton *whatsThisButton = new QPushButton(mainGroup);
    whatsThisButton->setText(tr("What's This?"));
    whatsThisButton->setCheckable(true);

    QHBoxLayout *viewLayout = new QHBoxLayout(this);
    //viewLayout->addWidget(view);
    viewLayout->addWidget(mainGroup);

    QVBoxLayout *rotateGroupLayout = new QVBoxLayout(rotateGroup);
    rotateGroupLayout->addWidget(rotateSlider);

    QVBoxLayout *scaleGroupLayout = new QVBoxLayout(scaleGroup);
    scaleGroupLayout->addWidget(scaleSlider);

    QVBoxLayout *shearGroupLayout = new QVBoxLayout(shearGroup);
    shearGroupLayout->addWidget(shearSlider);

    QVBoxLayout *typeGroupLayout = new QVBoxLayout(typeGroup);
    typeGroupLayout->addWidget(vectorType);
    typeGroupLayout->addWidget(pixmapType);
    typeGroupLayout->addWidget(textType);
    typeGroupLayout->addSpacing(4);
    typeGroupLayout->addWidget(textEditor);

    QVBoxLayout *mainGroupLayout = new QVBoxLayout(mainGroup);
    mainGroupLayout->addWidget(rotateGroup);
    mainGroupLayout->addWidget(scaleGroup);
    mainGroupLayout->addWidget(shearGroup);
    mainGroupLayout->addWidget(typeGroup);
    mainGroupLayout->addStretch(1);
    mainGroupLayout->addWidget(resetButton);
    mainGroupLayout->addWidget(animateButton);
    mainGroupLayout->addWidget(showSourceButton);
#ifdef QT_OPENGL_SUPPORT
    mainGroupLayout->addWidget(enableOpenGLButton);
#endif
    mainGroupLayout->addWidget(whatsThisButton);

    /*connect(rotateSlider, SIGNAL(valueChanged(int)), view, SLOT(changeRotation(int)));
    connect(shearSlider, SIGNAL(valueChanged(int)), view, SLOT(changeShear(int)));
    connect(scaleSlider, SIGNAL(valueChanged(int)), view, SLOT(changeScale(int)));

    connect(vectorType, SIGNAL(clicked()), view, SLOT(setVectorType()));
    connect(pixmapType, SIGNAL(clicked()), view, SLOT(setPixmapType()));
    connect(textType, SIGNAL(clicked()), view, SLOT(setTextType()));
    connect(textType, SIGNAL(toggled(bool)), textEditor, SLOT(setEnabled(bool)));
    connect(textEditor, SIGNAL(textChanged(QString)), view, SLOT(setText(QString)));

    connect(view, SIGNAL(rotationChanged(int)), rotateSlider, SLOT(setValue(int)));
    connect(view, SIGNAL(scaleChanged(int)), scaleSlider, SLOT(setValue(int)));
    connect(view, SIGNAL(shearChanged(int)), shearSlider, SLOT(setValue(int)));*/

    /*connect(resetButton, SIGNAL(clicked()), view, SLOT(reset()));
    connect(animateButton, SIGNAL(clicked(bool)), view, SLOT(setAnimation(bool)));
    connect(whatsThisButton, SIGNAL(clicked(bool)), view, SLOT(setDescriptionEnabled(bool)));
    connect(whatsThisButton, SIGNAL(clicked(bool)), view->hoverPoints(), SLOT(setDisabled(bool)));
    connect(view, SIGNAL(descriptionEnabledChanged(bool)), view->hoverPoints(), SLOT(setDisabled(bool)));
    connect(view, SIGNAL(descriptionEnabledChanged(bool)), whatsThisButton, SLOT(setChecked(bool)));
    connect(showSourceButton, SIGNAL(clicked()), view, SLOT(showSource()));
#ifdef QT_OPENGL_SUPPORT
    connect(enableOpenGLButton, SIGNAL(clicked(bool)), view, SLOT(enableOpenGL(bool)));
#endif
    view->loadSourceFile(":res/ffans/control.cpp");
    view->loadDescription(":res/ffans/control.html");*/

    // defaults
    //view->reset();
    vectorType->setChecked(true);
    textEditor->setText("Qt Affine Transformation Example");
    textEditor->setEnabled(false);

    //animateButton->animateClick();
}
