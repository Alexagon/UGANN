#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    qDebug () << "It's okey";
    mw2 = new MainWindow2;
    qsrand((uint)QTime::currentTime().msec());

    ui->pushButton_2->setEnabled(false);

    ui->lineEdit_2->setText(QString::number(GA_in));
    ui->lineEdit_3->setText(QString::number(GA_out));
    ui->lineEdit_4->setText(QString::number(GA_layers));
    ui->lineEdit_5->setText(QString::number(GA_lel));
    ui->lineEdit_6->setText(QString::number(GA_pop));
    ui->lineEdit_7->setText(QString::number(GA_bw));
    ui->lineEdit_8->setText(QString::number(GA_crosschance));
    ui->lineEdit_9->setText(QString::number(GA_mutchance));
    ui->lineEdit_11->setText(QString::number(GA_allowed_generations));


    //==========================================================================================
    //QFile file("wpbc(new).data");
    //QFile file("Breast Cancer Wisconsin (Diagnostic) Data Set (for reading).txt");
    QFile file("breast-cancer-wisconsin (for read).txt");

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        QMessageBox::information(this, "error", file.errorString());

    QStringList lines, line;
    int counter = 1;

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString l = in.readLine();
        lines.append(l);
        counter++;
    }
    file.close();
    counter--;
    qDebug () << "we cool?" << counter;

    QVector <int> remove_this;
    remove_this.append(1);

    int pp = remove_this.size();

    // ======================================================= ПРИЗНАКИ НА ВЫБРОС ЗА CART
    //remove_this.append({5, 6, 11, 34, 35});// --------------- == 1 == (1) wpbc(new).data
    //remove_this.append({8, 9, 10, 13});// ------------------- == 2 == (1) wpbc(new).data
    //remove_this.append({5, 6, 8, 10, 11, 34, 35});// -------- = 1.5 = (1) wpbc(new).data
    //remove_this.append({3, 12, 11, 8, 9});// ---------------- == 1 == (2) Breast Cancer Wisconsin (Diagnostic) Data Set (for reading).txt
    remove_this.append({7, 9, 10, 11});// ------------------- == 1 == (3) breast-cancer-wisconsin (for read).txt

    pp = remove_this.size() - pp;
    ui->lineEdit_2->setText(QString::number(GA_in-pp));

    std::sort(remove_this.begin(), remove_this.end());
    GA_data.clear();
    for (int i = 0; i < lines.size(); i++)
    {
        line = lines[i].split(',');
        for (int j = remove_this.size()-1; j >= 0; j--)
            line.removeAt(remove_this[j]-1); // ============= ОТБОР ПРИЗНАКОВ ЗА CART
        GA_data.append(QVector<double>());
        for (int j = 0 ; j < line.size(); j++)
            GA_data.last().append(line[j].toDouble());
    }
    qDebug () << "Data is opened\nAnd its size:" << GA_data.size() << GA_data.first().size();
    double cur_class_coef, class_coef = 0;
    int data_num, ti;
    for (int i = 0; i < GA_data.size(); i++)
        class_coef += GA_data[i][0];
    class_coef /= GA_data.size();
    cur_class_coef = class_coef;
    data_num = int(double(1-GA_split_data_coef)*GA_data.size())+1;
    for (int i = 0; i < data_num; i++)
    {
        if (cur_class_coef == class_coef)
            ti = qrand () % GA_data.size();
        if (cur_class_coef < class_coef)
        {
            do {
                ti = qrand () % GA_data.size();
            } while (GA_data[ti][0] != 1);
        }
        if (cur_class_coef > class_coef)
        {
            do {
                ti = qrand () % GA_data.size();
            } while (GA_data[ti][0] != 0);
        }
        GA_exp_data.append(GA_data.takeAt(ti));
        cur_class_coef = 0;
        for (int j = 0; j < GA_exp_data.size(); j++)
            cur_class_coef += GA_exp_data[j][0];
        cur_class_coef /= GA_exp_data.size();
    }
    qDebug () << "Number" << GA_split_data_coef << GA_data.size() << GA_exp_data.size() << GA_data.size()+GA_exp_data.size();
    qDebug () << "Coef" << class_coef << cur_class_coef;
    //==========================================================================================
}

MainWindow::~MainWindow()
{
    delete ui;
}


QVector <int> MainWindow::decoder(QString chr)
{
    QVector <int> v;
    int a;
    int helper;
    for (int i = 0; i < GA_links; i++)
    {
        a = 0;
        for (int j = 1; j < GA_bw; j++)
        {
            helper = GA_bw - j - 1;
            a += chr.mid(i*GA_bw+j, 1).toInt() * qPow(2, helper);
        }
        a = (bool(chr.mid(i*GA_bw, 1).toInt()))?a:-a;
        v.append(a);
    }
    return v;
}

double MainWindow::calc_errors(QVector<int> v, bool is_best)
{
    QVector <QVector <double> > in;
    QVector <double> valid, clas;
    QVector <double> r;
    double x;

    if (!is_best)
        in = GA_data;
    else
        in = GA_exp_data;

    QVector <double> last, input, output;
    QVector <double> truth;
    int tmp, pointer, errors = 0;
    double avg_error;
    r.append(GA_in);
    for (int i = 0; i < GA_layers; i++)
        r.append(GA_lel);
    r.append(GA_out);
    for (int i = 0; i < in.size(); i++)
        truth.append(in[i].takeAt(0));
    for (int i = 0; i < in.size(); i++) // Входы
    {
        QApplication::processEvents();
        pointer = 0;
        last.clear();
        for (int j = 0; j < r.size()-1; j++) // Этапы перцептрона
        {
            input.clear();
            input = (j==0)?in[i]:output;
            output.clear();
            for (int k = 0; k < r[j+1]; k++) // Output
            {
                tmp = 0;
                for (int l = 0; l < r[j]; l++) // Input
                    tmp += (input[l] * v[pointer + l*r[j+1]]); // Умножение на веса
                tmp += v[pointer+r[j]*r[j+1]+k]; // Добавляем сдвиг
                if (j+1 != r.size()-1)
                    output.append(1 / (1 + qExp(double(-tmp) * 2)));
                else
                {
                    last.append((tmp>0)?1:0);
                    if (is_best)
                    {
                        valid.append(double(tmp));
                        clas.append(last.first());
                    }
                }
            }
            pointer += ((r[j]+1)*r[j+1]);
        }
        input.clear();
        output.clear();
        if ((is_best) && (x != in[i][0]))
            x = in[i][0];
        errors += qPow((truth[i] - last.first()), 2); // --------------------- ТОЛЬКО ДЛЯ ОДНОГО ВЫХОДА
    }
    avg_error = double(errors) / double(in.size());
    double valid_min, valid_max;
    if (is_best)
    {
        valid_min = valid.first();
        valid_max = valid.first();
        for (int i = 0; i < valid.size(); i++)
        {
            if (valid[i] > valid_max)
                valid_max = valid[i];
            if (valid[i] < valid_min)
                valid_min = valid[i];
        }
        for (int i = 0; i < valid.size(); i++)
        {
            if (valid[i] > 0)
                valid[i] = valid[i]/qAbs(valid_max);
            else
                valid[i] = qAbs(valid[i])/qAbs(valid_min);
            //qDebug () << i+1 << clas[i] << valid[i];
            std::cout << i+1 << "	" << clas[i] << "	" << valid[i] << std::endl;
        }
    }
    return avg_error;
}

void MainWindow::crossover(QString a, QString b)
{
    QString kid1, kid2;
    int p = qrand () % (GA_links - 1);
    kid1 = a.mid(0, (p+1)*GA_bw) + b.mid((p+1)*GA_bw, (GA_links-p-1)*GA_bw);
    kid2 = b.mid(0, (p+1)*GA_bw) + a.mid((p+1)*GA_bw, (GA_links-p-1)*GA_bw);
    GA_ppl.append(chromosome(kid1));
    GA_ppl.append(chromosome(kid2));
}

void MainWindow::mutation(QString m)
{
    int p = qrand () % m.size();
    m[p] = (bool(m[p].digitValue()) == false)?'1':'0';
    GA_ppl.append(chromosome(m));
}

void MainWindow::create_gen()
{
    GA_ppl.clear();
    for (int i = 0; i < GA_pop; i++)
    {
        GA_ppl.append(chromosome(""));
        for (int j = 0; j < GA_bits; j++)
            GA_ppl[i].seq.append(QString::number(qrand () % 2));
    }
    generation++;
    ui->label_2->setText("GENERATION: " + QString::number(generation));
}

void MainWindow::evaluation()
{
    int done;
    bool ibc = false;
    ui->progressBar->setValue(0);
    for (int i = 0; (i < GA_ppl.size() && !stop); i++)
    {
        GA_ppl[i].error = calc_errors(decoder(GA_ppl[i].seq), false);
        if (GA_ppl[i].error < best.error)
        {
            best.seq = GA_ppl[i].seq;
            best.error = GA_ppl[i].error;
            ibc = true;
        }
        done = (i+1)*100/GA_ppl.size();
        ui->label->setText(QString::number(done) + "% Done (" + QString::number(i+1) + "/" + QString::number(GA_ppl.size()) + ")");
        ui->progressBar->setValue(done);
    }
    if (ibc)
    {
        ui->lineEdit->setText("BEST: " + QString::number(best.error) + " (G-" + QString::number(generation) + ")");
        QApplication::alert(this, 0);
    }
}

void MainWindow::selection()
{
    if (!stop)
    {
        std::sort(GA_ppl.begin(), GA_ppl.end());
        if (GA_ppl.size() > GA_pop)
            GA_ppl.erase(GA_ppl.begin()+GA_pop, GA_ppl.end());
        for (int i = 3; i < GA_ppl.size(); i++)
            GA_ppl.removeAt(i);
    }
}

void MainWindow::reproduction()
{
    int size = GA_ppl.size();
    int p1, p2, m;

    for (int i = 0; i < (size*GA_crosschance/100); i++)
    {
        p1 = qrand () % size;
        do
            p2 = qrand () % size;
        while (p1 == p2);
        crossover(GA_ppl[p1].seq, GA_ppl[p2].seq);
    }

    for (int i = 0; i < (size*GA_mutchance/100); i++)
    {
        m = qrand () % size;
        mutation(GA_ppl[m].seq);
    }

    do
        p2 = qrand () % size;
    while (p2 == 0);
    crossover(GA_ppl[0].seq, GA_ppl[p2].seq);
    mutation(GA_ppl[0].seq);

    GA_ppl.erase(GA_ppl.begin()+1, GA_ppl.begin()+size);
    generation++;
    ui->label_2->setText("GENERATION: " + QString::number(generation));
}

void MainWindow::send_results(QString res)
{
    QVector <int> v;
    v = decoder(best.seq);
    best.error = calc_errors(v, true);
    mw2->set_results(best, v, GA_bw, res, generation, GA_resolution);
    mw2->show_results();
    mw2->show();
    QApplication::alert(mw2, 0);
    this->close();
}

void MainWindow::on_pushButton_clicked()
{
    bool ready = true;
    if ((ui->lineEdit_2->text().isEmpty()) || (ui->lineEdit_3->text().isEmpty()) || (ui->lineEdit_4->text().isEmpty())
            || (ui->lineEdit_5->text().isEmpty()) || (ui->lineEdit_6->text().isEmpty()) || (ui->lineEdit_7->text().isEmpty())
            || (ui->lineEdit_8->text().isEmpty()) || (ui->lineEdit_9->text().isEmpty())
            || (ui->lineEdit_11->text().isEmpty()))
        ready = false;
    if (ready)
    {
        GA_in = ui->lineEdit_2->text().toInt();
        GA_out = ui->lineEdit_3->text().toInt();
        GA_layers = ui->lineEdit_4->text().toInt();
        GA_lel = ui->lineEdit_5->text().toInt();
        GA_pop = ui->lineEdit_6->text().toInt();
        GA_bw = ui->lineEdit_7->text().toInt();
        GA_crosschance = ui->lineEdit_8->text().toInt();
        GA_mutchance = ui->lineEdit_9->text().toInt();
        GA_target = 0;
        GA_allowed_generations = ui->lineEdit_11->text().toInt();
        GA_links = GA_in*GA_lel + (GA_layers-1)*GA_lel*GA_lel + GA_layers*GA_lel + GA_out*(GA_lel+1);
        GA_bits = GA_links * GA_bw;
        GA_cb = 0;

        ui->lineEdit_2->setEnabled(false);
        ui->lineEdit_3->setEnabled(false);
        ui->lineEdit_4->setEnabled(false);
        ui->lineEdit_5->setEnabled(false);
        ui->lineEdit_6->setEnabled(false);
        ui->lineEdit_7->setEnabled(false);
        ui->lineEdit_8->setEnabled(false);
        ui->lineEdit_9->setEnabled(false);
        ui->lineEdit_11->setEnabled(false);

        ui->pushButton_2->setEnabled(true);

        create_gen();
        evaluation();
        selection();
        while ((generation < GA_allowed_generations) && (best.error > GA_target) && !stop)
        {
            reproduction();
            evaluation();
            selection();
        }
        if (generation >= GA_allowed_generations)
            send_results("The algorithm is finished (reached limit of generations)");
        if (best.error <= GA_target)
            send_results("The algorithm is finished (found the optimal solution)");
    }
    else
        QMessageBox::warning(this, "Whoops", "You need to fill all the fields to start");
}

void MainWindow::on_pushButton_2_clicked()
{
    QMessageBox::StandardButton ask;
    ask = QMessageBox::question(this, "Wait", "Do you want to see the results, "
                                              "with the current best solution?",
                                QMessageBox::No|QMessageBox::Yes|QMessageBox::Cancel);
    if (ask == QMessageBox::Yes)
    {
        stop = true;
        send_results("The algorithm was canceled");
    }
    else if (ask == QMessageBox::No)
        exit(EXIT_SUCCESS);
}
