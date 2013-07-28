#ifndef MODELGRAPH_H
#define MODELGRAPH_H

#include <QObject>
#include <QColor>
//#include <QGLWidget>

class Patch;
struct Geometry;

//! [0]
class ModelGraph : public QObject
{
public:
    explicit ModelGraph(QObject *parent, int d = 64, qreal s = 1.0);
    ~ModelGraph();
    void setColor(QColor c);
    void draw() const;
private:
    void buildGeometry(int d, qreal s);

    QList<Patch *> parts;
    Geometry *geom;

    QColor cartesianSystemColor;
};
//! [0]

#endif // MODELGRAPH_H
