SOURCES += \
    glwidget.cpp \
    main.cpp \
    modelwindow.cpp \
    modeldb.cpp \
    modelsimcore.cpp \
    modeldataanalysis.cpp \
    settingswindow.cpp \
    settingsxml.cpp \
    modeldataanalysiswindow.cpp \
    modelparameterswindow.cpp \
    modelparameters.cpp \
    controlwidget.cpp \
    modelgraphics.cpp
HEADERS += \
    glwidget.h \
    main.h \
    modelwindow.h \
    controlwidget.h \
    modeldb.h \
    modelsimcore.h \
    modeldataanalysis.h \
    settingswindow.h \
    settingsxml.h \
    modeldataanalysiswindow.h \
    modelparameterswindow.h \
    modelparameters.h \
    modelgraphics.h

qtHaveModule(opengl) {
	DEFINES += QT_OPENGL_SUPPORT
	QT += opengl
}
QT += widgets

#SHARED_FOLDER = ../shared

#include($$SHARED_FOLDER/shared.pri)

#RESOURCES += \
#    ffans.qrc

# install
#target.path = $$[QT_INSTALL_EXAMPLES]/widgets/painting/affine
#INSTALLS += target

#wince*: {
#    DEPLOYMENT_PLUGIN += qjpeg
#}
