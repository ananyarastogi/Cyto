// 21-12-2017

#ifndef DISTRIBWIDGET_H
#define DISTRIBWIDGET_H

#include <QWidget>
#include <QListWidget>
#include "distribution.h"
#include "random.h"
//#include "viewdistrib.h"
struct viewDistrib; // avoids circular reference ?


/// This class maintains the currentLaw 'Law' according to the options.

namespace Ui {
class DistribWidget;
}

class DistribWidget : public QWidget
{
    Q_OBJECT


public:
    // be careful to have same order in the combo and enum
    explicit DistribWidget(QWidget *parent = 0);
    ~DistribWidget();
    void Reset();
    void enableCombo(int index, bool enable);
    int nbCombos();

    QString currentFile;
    QListWidget* listCombo;
    Law* getLaw();
    void setLaw(Law & toCopy);

    // only update the spin boxes, doesn't change the Law
    void setMu1(double _mu1);
    void setMu2(double _mu2);
    void setSigma1(double _Sigma1);
    void setSigma2(double _Sigma2);
    void setWeight(double _weight);
    void setDistrib(int typeDistrib);

public slots:
    void ComboChanged(bool allowLoading = true);
    void LoadData(QString fileToOpen);
    void update();
    void ShowPlot();
    void restoreParent(viewDistrib *toBeDisconnected);

private:
    QWidget* initialParent;
    Law* currentLaw;
    Ui::DistribWidget *ui;
};

#endif // DISTRIBWIDGET_H
