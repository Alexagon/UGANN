#ifndef MAINWINDOW2_H
#define MAINWINDOW2_H

#include <QMainWindow>
#include <QDebug>
#include <QString>
#include <QVector>

#include "chromosome.h"

namespace Ui {
class MainWindow2;
}

class MainWindow2 : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow2(QWidget *parent = nullptr);
    ~MainWindow2();

    void set_results(chromosome chr, QVector <int> vec, int b_weight,
                     QString res, int gen, int resol);
    void show_results();

private:
    Ui::MainWindow2 *ui;

    chromosome best = chromosome ("", 1);
    QVector <int> v;
    int bw;
    QString result;
    int generations;
    int resolution;

    bool vis = false;

};

#endif // MAINWINDOW2_H
