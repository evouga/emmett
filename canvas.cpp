#include "canvas.h"
#include "mainwindow.h"

#include <QTimer>

Canvas::Canvas(QWidget *parent) :
    QGLWidget(parent), mw_(NULL), tickcount_(0)
{
    timer_ = new QTimer(this);
    connect(timer_, SIGNAL(timeout()), this, SLOT(tick()));
    timer_->start(33);
}

Canvas::~Canvas()
{
    delete timer_;
}

void Canvas::setMainWindow(MainWindow *mw)
{
    mw_ = mw;
    centerCamera();
}

void Canvas::centerCamera()
{
    assert(mw_);
    mw_->getSceneBounds(center_, radius_);
}

void Canvas::initializeGL()
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glDisable(GL_DEPTH_TEST);
}

void Canvas::resizeGL(int w, int h)
{
    glViewport(0,0,w,h);
}

void Canvas::paintGL()
{
    assert(mw_);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(center_[0]-radius_, center_[0]+radius_, center_[1]-radius_, center_[1]+radius_, 0, 1);
    glMatrixMode(GL_MODELVIEW);

    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear (GL_COLOR_BUFFER_BIT);
    glColor3f (0.0, 0.0, 0.0);

    mw_->renderScene(tickcount_);
}

void Canvas::tick()
{
    tickcount_++;
    updateGL();
}
