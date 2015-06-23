#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <Eigen/Core>

class Morphing;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void renderScene(int tick);
    void getSceneBounds(Eigen::Vector2d &center, double &radius);
    
private slots:
    void on_actionLoad_Shapes_triggered();

private:
    bool parseGeometryFile(std::string filename);
    void renderPolygon(const Eigen::VectorXd &polygon);

    Ui::MainWindow *ui;   

    Eigen::VectorXd polygon1_;
    Eigen::VectorXd polygon2_;
    Morphing *morph_;
};

#endif // MAINWINDOW_H
