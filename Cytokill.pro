

HEADERS += \
    main.h \
    random.h \
    cell.h \
    interact.h \
    distribution.h \
    hypothesis.h \
    outsidefunc.h \
    plot3d.h \
    starter.h \
    output.h \
    dataStorage.h \
    observer.h \
    trackball.h \
    distribwidget.h \
    viewdistrib.h \
    grapheCustom.h \
    QCustomPlot/qcustomplot.h \
    parameters.h


SOURCES += \
    main.cpp \
    random.cpp \
    cell.cpp \
    interact.cpp \
    distribution.cpp \
    hypothesis.cpp \
    outsidefunc.cpp \
    plot3d.cpp \
    starter.cpp \
    output.cpp \
    dataStorage.cpp \
    observer.cpp \
    trackball.cpp \
    distribwidget.cpp \
    viewdistrib.cpp \
    grapheCustom.cpp \
    QCustomPlot/qcustomplot.cpp \
    parameters.cpp

QMAKE_CXXFLAGS += -std=c++11

#this is required for qt to create the ui_....h automatically
QT += widgets printsupport svg core gui

FORMS += \
    starter.ui \
    distribwidget.ui \
    viewdistrib.ui

unix: LIBS += -lglut -lGLU -lGL -lSOIL

win32: LIBS += -L$$PWD/freeglut/lib/ -lfreeglut -lopengl32
INCLUDEPATH += $$PWD/freeglut/include
#DEPENDPATH += $$PWD/freeglut/lib

win32: LIBS += -L$$PWD/soil/ -lSOIL
INCLUDEPATH += $$PWD/soil
DEPENDPATH += $$PWD/soil

win32: LIBS += -L$$PWD/glStatic/ -lglu32
INCLUDEPATH += $$PWD/glStatic
DEPENDPATH += $$PWD/glStatic

win32: LIBS += -L$$PWD/glStatic/ -lopengl32
INCLUDEPATH += $$PWD/glStatic
DEPENDPATH += $$PWD/glStatic
