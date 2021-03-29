#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include <QBitArray>
#include <QTime>
#include <QStringList>
#include <QVector>
#include <QList>
#include <QtMath>
#include <QMessageBox>
#include <QTextStream>
#include <QFile>
#include <iostream>
#include <algorithm>

#include "chromosome.h"
#include <mainwindow2.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    QVector <int> decoder (QString chr);
    double calc_errors (QVector <int>, bool is_best);
    void crossover (QString a, QString b);
    void mutation (QString m);
    void create_gen ();
    void evaluation ();
    void selection ();
    void reproduction ();
    void send_results(QString res);

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();


private:
    Ui::MainWindow *ui;

    int GA_in = 9;       // Входы
    int GA_out = 1;      // Выходы

    int GA_layers = 3;   // Скрытые слои
    int GA_lel = 20;      // Элементов в одном слое

    int GA_pop = 1000;     // Популяция
    int GA_bw = 7;       // Колличество битов в весах

    int GA_crosschance = 100;  // Шанс кроссовера
    int GA_mutchance = 10;     // Шанс мутации

    double GA_target = 0;                // Целевая ошибка ГА
    int GA_allowed_generations = 1000;      // Количество позволенных поколений

    int GA_resolution = 10; // [0 .. 100] Размерность плоскости
    int GA_step = 1; // Шаг по плоскости для точек обучения

    double GA_split_data_coef = 0.75;

    int GA_links = GA_in*GA_lel + (GA_layers-1)*GA_lel*GA_lel + GA_layers*GA_lel + GA_out*(GA_lel+1); // Колличество связей
    int GA_bits = GA_links * GA_bw;         // Колличество битов в НС
    chromosome best = chromosome("", 1);
    int generation = 0;

    QList <chromosome> GA_ppl;

    bool stop = false;
    int GA_cb;
    QVector <QVector<double>> GA_data, GA_exp_data;


    MainWindow2 *mw2;
};

#endif // MAINWINDOW_H

