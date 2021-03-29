#include "mainwindow2.h"
#include "ui_mainwindow2.h"

MainWindow2::MainWindow2(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow2)
{
    ui->setupUi(this);
}

MainWindow2::~MainWindow2()
{
    delete ui;
}

void MainWindow2::set_results(chromosome chr, QVector <int> vec, int b_weight,
                              QString res, int gen, int resol)
{
    best.seq = chr.seq;
    best.error = chr.error;
    v = vec;
    bw = b_weight;
    result = res;
    generations = gen;
    resolution = resol;
}

void MainWindow2::show_results()
{
    ui->lineEdit->setText(result);
    ui->lineEdit_2->setText(QString::number(generations));
    ui->lineEdit_3->setText(QString::number(best.error));
    QString t = "";
    for (int i = 0; i < v.size(); i++)
        t.append(" " + best.seq.mid(i*bw, 1) + " " + best.seq.mid(i*bw+1, bw-1) +
                 ((v[i]>=0)?"               ":"              ") +
                 QString::number(v[i]) + "\n");
    t.chop(1);
    ui->textEdit->setText(t);
}
