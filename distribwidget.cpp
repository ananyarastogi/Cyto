#include "distribwidget.h"
#include "ui_distribwidget.h"
#include <QStringList>
#include <QFileDialog>
#include "viewdistrib.h"

DistribWidget::DistribWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DistribWidget)
{
    ui->setupUi(this);

    initialParent = parent;
    listCombo = NULL;

    QStringList L = {QString("Fixed"), QString("Normal"), QString("LogNormal"), QString("From Data"), QString("BiModal(Normal)"), QString("Exponential")};
    if(L.size() != NBDIstribs) cerr << "You NERD, fix this ! Distrib enum != combo size" << endl;

    // (normal way to make the combo box) ui->comboBoxTypeDistrib->addItems(L);
    // more tricks here: http://stackoverflow.com/questions/11439773/disable-item-in-qt-combobox
    // item->setFlags(disable ? item->flags() & ~(Qt::ItemIsSelectable|Qt::ItemIsEnabled)
    // : Qt::ItemIsSelectable|Qt::ItemIsEnabled));
    // => The way to make a combo box with hideable items :
    listCombo = new QListWidget();//ui->comboBoxTypeDistrib);    // the combo box for the experiments is made using ListWidgets (instead of usual combo), because they can be enabled / disabled
    //listCombo->hide();
    listCombo->addItems(L);
    ui->comboBoxTypeDistrib->setModel(listCombo->model());
    //ui->comboBoxTypeDistrib->clear();
   /* for(int i = 0; i < (int) listCombo->count(); ++i){
        QListWidgetItem *item = listCombo->item(i);
        item->setFlags(item->flags()); // & Qt::ItemIsEnabled);
    }*/
    //QObject::connect(ui->comboBoxTypeDistrib, SIGNAL(activated(int)), this, SLOT(ComboChanged()));
    //ui->comboBoxTypeDistrib->setCurrentIndex(0); // wonder what happens when item 0 is not enabled, later ...



    currentLaw = new Law();
    Reset();
    QObject::connect(ui->comboBoxTypeDistrib, SIGNAL(activated(int)), this, SLOT(ComboChanged()));
    QObject::connect(ui->doubleSpinBoxMu1, SIGNAL(valueChanged(double)), this, SLOT(update()));
    QObject::connect(ui->doubleSpinBoxSigma1, SIGNAL(valueChanged(double)), this, SLOT(update()));
    QObject::connect(ui->doubleSpinBoxWeight, SIGNAL(valueChanged(double)), this, SLOT(update()));
    QObject::connect(ui->doubleSpinBoxMu2, SIGNAL(valueChanged(double)), this, SLOT(update()));
    QObject::connect(ui->doubleSpinBoxSigma2, SIGNAL(valueChanged(double)), this, SLOT(update()));
    QObject::connect(ui->pushButton, SIGNAL(released()), this, SLOT(ShowPlot()));
    ComboChanged();
}

void DistribWidget::setLaw(Law & toCopy){
    currentLaw = new Law(toCopy); // better to copy from a simulation
    setDistrib(toCopy.type);
    setMu1(toCopy.mu1);
    setMu2(toCopy.mu2);
    setSigma1(toCopy.sigma1);
    setSigma2(toCopy.sigma2);
    setWeight(toCopy.weight);
    currentFile = QString(toCopy.sourceDataFile.c_str());
    ComboChanged(false);
}


void DistribWidget::setDistrib(int typeDistrib){
    if((typeDistrib < 0) || (typeDistrib >= ui->comboBoxTypeDistrib->count())){
        cerr << "ERR: DistribWidget::setDistrib(" << typeDistrib << "), index out of bounds " << endl;
    }
    ui->comboBoxTypeDistrib->setCurrentIndex(typeDistrib);
}

void DistribWidget::restoreParent(viewDistrib* toBeDisconnected){
    QObject::disconnect(ui->comboBoxTypeDistrib, SIGNAL(activated(int)), toBeDisconnected, SLOT(replot()));
    QObject::disconnect(ui->doubleSpinBoxMu1, SIGNAL(valueChanged(double)), toBeDisconnected, SLOT(replot()));
    QObject::disconnect(ui->doubleSpinBoxSigma1, SIGNAL(valueChanged(double)), toBeDisconnected, SLOT(replot()));
    QObject::disconnect(ui->doubleSpinBoxWeight, SIGNAL(valueChanged(double)), toBeDisconnected, SLOT(replot()));
    QObject::disconnect(ui->doubleSpinBoxMu2, SIGNAL(valueChanged(double)), toBeDisconnected, SLOT(replot()));
    QObject::disconnect(ui->doubleSpinBoxSigma2, SIGNAL(valueChanged(double)), toBeDisconnected, SLOT(replot()));

    this->setParent(initialParent);
    this->show();
    ui->pushButton->show();
}

void DistribWidget::ShowPlot(){
    ui->pushButton->hide();
    //QWidget* previousParent = (QWidget*) this->parent();
    viewDistrib* a = new viewDistrib(this);//doesn't work with non-pointers, weird.

    QObject::connect(ui->comboBoxTypeDistrib, SIGNAL(activated(int)), a, SLOT(replot()));
    QObject::connect(ui->doubleSpinBoxMu1, SIGNAL(valueChanged(double)), a, SLOT(replot()));
    QObject::connect(ui->doubleSpinBoxSigma1, SIGNAL(valueChanged(double)), a, SLOT(replot()));
    QObject::connect(ui->doubleSpinBoxWeight, SIGNAL(valueChanged(double)), a, SLOT(replot()));
    QObject::connect(ui->doubleSpinBoxMu2, SIGNAL(valueChanged(double)), a, SLOT(replot()));
    QObject::connect(ui->doubleSpinBoxSigma2, SIGNAL(valueChanged(double)), a, SLOT(replot()));

    a->show();
    //this->setParent(previousParent);
}

void DistribWidget::setMu1(double _mu1){
    ui->doubleSpinBoxMu1->setValue((std::isnan(_mu1)) ? 0 : _mu1);
}

void DistribWidget::setMu2(double _mu2){
    ui->doubleSpinBoxMu2->setValue((std::isnan(_mu2)) ? 0 : _mu2);
}

void DistribWidget::setSigma1(double _Sigma1){
    ui->doubleSpinBoxSigma1->setValue((std::isnan(_Sigma1)) ? 0 : _Sigma1);
}

void DistribWidget::setSigma2(double _Sigma2){
    ui->doubleSpinBoxSigma2->setValue((std::isnan(_Sigma2)) ? 0 : _Sigma2);
}

void DistribWidget::setWeight(double _weight){
    ui->doubleSpinBoxWeight->setValue((std::isnan(_weight)) ? 0 : _weight);
}


void DistribWidget::Reset(){
    currentLaw->set(Fixed, 0.0);
    currentFile = QString("");
    ui->doubleSpinBoxMu1->setValue(10);
    ui->doubleSpinBoxSigma1->setValue(2);
    ui->doubleSpinBoxMu2->setValue(30);
    ui->doubleSpinBoxSigma2->setValue(1);
    ui->doubleSpinBoxWeight->setValue(0.8);
}

DistribWidget::~DistribWidget()
{
    cerr << "Destroy ???" << endl;
    delete ui;
}

void DistribWidget::ComboChanged(bool allowLoading){
    switch(ui->comboBoxTypeDistrib->currentIndex()){
    case Fixed:{
        ui->doubleSpinBoxMu1->show();
        ui->doubleSpinBoxSigma1->hide();
        ui->doubleSpinBoxMu2->hide();
        ui->doubleSpinBoxSigma2->hide();
        ui->doubleSpinBoxWeight->hide();
        ui->labelMean->show();
        ui->labelStdDev->hide();
        ui->labelMean2->hide();
        ui->labelStdDev2->hide();
        ui->labelWeight->hide();
        ui->pushButton->hide();
        break;}
    case Normal:{
        ui->doubleSpinBoxMu1->show();
        ui->doubleSpinBoxSigma1->show();
        ui->doubleSpinBoxMu2->hide();
        ui->doubleSpinBoxSigma2->hide();
        ui->doubleSpinBoxWeight->hide();
        ui->labelMean->show();
        ui->labelStdDev->show();
        ui->labelMean2->hide();
        ui->labelStdDev2->hide();
        ui->labelWeight->hide();
        ui->pushButton->show();
        break;}
    case LogNormal:{
        ui->doubleSpinBoxMu1->show();
        ui->doubleSpinBoxSigma1->show();
        ui->doubleSpinBoxMu2->hide();
        ui->doubleSpinBoxSigma2->hide();
        ui->doubleSpinBoxWeight->hide();
        ui->labelMean->show();
        ui->labelStdDev->show();
        ui->labelMean2->hide();
        ui->labelStdDev2->hide();
        ui->labelWeight->hide();
        ui->pushButton->show();
        break;}
    case FromData:{
        ui->doubleSpinBoxMu1->hide();
        ui->doubleSpinBoxSigma1->hide();
        ui->doubleSpinBoxMu2->hide();
        ui->doubleSpinBoxSigma2->hide();
        ui->doubleSpinBoxWeight->hide();
        ui->labelMean->hide();
        ui->labelStdDev->hide();
        ui->labelMean2->hide();
        ui->labelStdDev2->hide();
        ui->labelWeight->hide();
        ui->pushButton->show();
        if(allowLoading){
            QString newFile = QFileDialog::getOpenFileName(this, currentFile);
            if(newFile.size() > 0) {
                currentFile = newFile;
                LoadData(currentFile);
            }
        }
        break;}
    case BiModal:{
        ui->doubleSpinBoxMu1->show();
        ui->doubleSpinBoxSigma1->show();
        ui->doubleSpinBoxMu2->show();
        ui->doubleSpinBoxSigma2->show();
        ui->doubleSpinBoxWeight->show();
        ui->labelMean->show();
        ui->labelStdDev->show();
        ui->labelMean2->show();
        ui->labelStdDev2->show();
        ui->labelWeight->show();
        ui->pushButton->show();
        break;}
    case Exponential:{
        ui->doubleSpinBoxMu1->show();
        ui->doubleSpinBoxSigma1->hide();
        ui->doubleSpinBoxMu2->hide();
        ui->doubleSpinBoxSigma2->hide();
        ui->doubleSpinBoxWeight->hide();
        ui->labelMean->show();
        ui->labelStdDev->hide();
        ui->labelMean2->hide();
        ui->labelStdDev2->hide();
        ui->labelWeight->hide();
        ui->pushButton->show();
        break;}
    default:{break;}
    }
    update();
}

void DistribWidget::update(){
    switch(ui->comboBoxTypeDistrib->currentIndex()){
    case Fixed:{currentLaw->set(Fixed, ui->doubleSpinBoxMu1->value()); break;}
    case Normal:{currentLaw->set(Normal, ui->doubleSpinBoxMu1->value(), ui->doubleSpinBoxSigma1->value()); break;}
    case LogNormal:{currentLaw->set(LogNormal, ui->doubleSpinBoxMu1->value(), ui->doubleSpinBoxSigma1->value()); break;}
    case FromData:{break;}
    case BiModal:{currentLaw->set(BiModal, ui->doubleSpinBoxMu1->value(), ui->doubleSpinBoxSigma1->value(), ui->doubleSpinBoxWeight->value(), ui->doubleSpinBoxMu2->value(), ui->doubleSpinBoxSigma2->value()); break;}
    case Exponential:{currentLaw->set(Exponential,ui->doubleSpinBoxMu1->value() ); break;}
    default:{break;}
    }
}

void DistribWidget::LoadData(QString fileToOpen){
    currentLaw->set(fileToOpen.toStdString());
}
Law *DistribWidget::getLaw(){
    return currentLaw;
}


void DistribWidget::enableCombo(int index, bool enable){
    if((index < 0) || (index >= (int) ui->comboBoxTypeDistrib->count())) {
        cerr << "ERR: DistribWidget::enableCombo(" << index << ",...), out of bounds (only " << ui->comboBoxTypeDistrib->count() << "combos" << endl; return;
    }
    //QObject::disconnect(ui->comboBoxTypeDistrib, SIGNAL(activated(int)), this, SLOT(ComboChanged()));

    QListWidgetItem *item = listCombo->item(index);
    // https://wiki.qt.io/QFlags_tutorial
    //item->setHidden(!enable);
    //else item->hide();
    //cout << ((enable) ? "Enable " : "Disable ") << "Distrib " << index << endl;
    bool disable = !enable;
    item->setFlags((disable) ? item->flags() & ~(Qt::ItemIsSelectable|Qt::ItemIsEnabled) : (Qt::ItemIsSelectable|Qt::ItemIsEnabled));
    //item->setHidden(disable); doesn't work
    // in case, change choice to first allowed item ?
    if(!(listCombo->item(ui->comboBoxTypeDistrib->currentIndex())->flags() & (Qt::ItemIsSelectable))){
        //cout << "Problem ! " << endl;
        for(int i = 0; i < (int) ui->comboBoxTypeDistrib->count(); ++i){
            if(listCombo->item(i)->flags() & (Qt::ItemIsSelectable)){
                //cout << i << " is good " << endl;
                ui->comboBoxTypeDistrib->setCurrentIndex(i);
                return;
            } //else cout << i << " is bad " << endl;
        }
    }

    //if(enable) {item->setFlags(item->flags() & ~Qt::ItemIsEnabled);}
    //else {item->setFlags(item->flags() & !Qt::ItemIsEnabled);}*/
    //QObject::connect(ui->comboBoxTypeDistrib, SIGNAL(activated(int)), this, SLOT(ComboChanged()));
}

int DistribWidget::nbCombos(){
    return ui->comboBoxTypeDistrib->count();
}

