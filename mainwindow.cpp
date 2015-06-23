#include "mainwindow.h"
#include <fstream>
#include <QMessageBox>
#include <QFileDialog>
#include "ui_mainwindow.h"
#include "morphing.h"
#include <iostream>

using namespace std;
using namespace Eigen;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow), morph_(NULL)
{
    ui->setupUi(this);
    ui->canvasWidget->setMainWindow(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::renderScene(int tick)
{
    glColor3f(1.0, 0.0, 0.0);
    renderPolygon(polygon1_);
    glColor3f(0.0, 1.0, 0.0);
    renderPolygon(polygon2_);

    if(morph_)
    {
        morph_->renderInterpolation(tick);
    }
}

void MainWindow::renderPolygon(const VectorXd &polygon)
{
    int npts = (int)polygon.size()/2;
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<npts; i++)
    {
        glVertex2d(polygon[2*i], polygon[2*i+1]);
    }
    if(npts > 0)
        glVertex2d(polygon[0], polygon[1]);
    glEnd();
}

void MainWindow::getSceneBounds(Eigen::Vector2d &center, double &radius)
{
    double minx = std::numeric_limits<double>::infinity();
    double miny = std::numeric_limits<double>::infinity();
    double maxx = -std::numeric_limits<double>::infinity();
    double maxy = -std::numeric_limits<double>::infinity();

    int pts1 = polygon1_.size()/2;
    int pts2 = polygon2_.size()/2;

    for(int i=0; i<pts1; i++)
    {
        minx = min(polygon1_[2*i], minx);
        miny = min(polygon1_[2*i+1], miny);
        maxx = max(polygon1_[2*i], maxx);
        maxy = max(polygon1_[2*i+1], maxy);
    }
    for(int i=0; i<pts2; i++)
    {
        minx = min(polygon2_[2*i], minx);
        miny = min(polygon2_[2*i+1], miny);
        maxx = max(polygon2_[2*i], maxx);
        maxy = max(polygon2_[2*i+1], maxy);
    }

    center[0] = 0.5*(minx+maxx);
    center[1] = 0.5*(miny+maxy);
    radius = max(maxx-minx, maxy-miny);
    radius *= 0.5 * 1.1;
}

void MainWindow::on_actionLoad_Shapes_triggered()
{
    QFileDialog savedialog(this, "Import Shape Geometry", ".", "Shape Files (*.geo)");
    savedialog.setFileMode(QFileDialog::AnyFile);
    savedialog.setDefaultSuffix("geo");
    savedialog.setViewMode(QFileDialog::List);
    if(savedialog.exec())
    {
        QStringList filenames = savedialog.selectedFiles();
        if(filenames.size() > 0)
        {
            QString filename = filenames[0];
            bool success = parseGeometryFile(filename.toStdString());
            if(!success)
            {
                QMessageBox errormsg;
                errormsg.setText(tr("Cannot open file ") + filename);
                errormsg.exec();
            }
            if(morph_)
            {
                delete morph_;
            }
            morph_ = new Morphing(MorphParams(), polygon1_, polygon2_);
            ui->canvasWidget->centerCamera();
        }
    }
}

bool MainWindow::parseGeometryFile(std::string filename)
{
    ifstream ifs(filename);
    if(!ifs)
        return false;

    int npts1;
    ifs >> npts1;
    if(!ifs)
        return false;

    VectorXd polygon1(npts1*2);
    for(int i=0; i<npts1; i++)
    {
        double x, y;
        ifs >> x >> y;
        polygon1[2*i] = x;
        polygon1[2*i+1] = y;
    }

    int npts2;
    ifs >> npts2;
    if(!ifs)
        return false;

    VectorXd polygon2(npts2*2);
    for(int i=0; i<npts2; i++)
    {
        double x, y;
        ifs >> x >> y;
        polygon2[2*i] = x;
        polygon2[2*i+1] = y;
    }

    if(!ifs)
        return false;

    polygon1_ = polygon1;
    polygon2_ = polygon2;
    return true;
}
