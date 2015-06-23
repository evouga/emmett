#ifndef CANVAS_H
#define CANVAS_H

#include <QGLWidget>
#include <Eigen/Core>

class MainWindow;

class Canvas : public QGLWidget
{
    Q_OBJECT
public:
    explicit Canvas(QWidget *parent = 0);
    virtual ~Canvas();

    void setMainWindow(MainWindow *mw);
    void centerCamera();

signals:
    
public slots:
    void tick();

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

private:

    MainWindow *mw_;
    Eigen::Vector2d center_;
    double radius_;
    QTimer *timer_;
    int tickcount_;
    
};

#endif // CANVAS_H
