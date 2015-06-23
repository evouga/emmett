#-------------------------------------------------
#
# Project created by QtCreator 2015-06-22T23:05:51
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

INCLUDEPATH += "C:/Program Files (x86)/Eigen_3.2.0"

TARGET = Morphing
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    morphing.cpp \
    canvas.cpp

HEADERS  += mainwindow.h \
    morphing.h \
    canvas.h

FORMS    += mainwindow.ui
