/****************************************************************************
** Meta object code from reading C++ file 'distribwidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.7.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "distribwidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'distribwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.7.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_DistribWidget_t {
    QByteArrayData data[11];
    char stringdata0[121];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_DistribWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_DistribWidget_t qt_meta_stringdata_DistribWidget = {
    {
QT_MOC_LITERAL(0, 0, 13), // "DistribWidget"
QT_MOC_LITERAL(1, 14, 12), // "ComboChanged"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 12), // "allowLoading"
QT_MOC_LITERAL(4, 41, 8), // "LoadData"
QT_MOC_LITERAL(5, 50, 10), // "fileToOpen"
QT_MOC_LITERAL(6, 61, 6), // "update"
QT_MOC_LITERAL(7, 68, 8), // "ShowPlot"
QT_MOC_LITERAL(8, 77, 13), // "restoreParent"
QT_MOC_LITERAL(9, 91, 12), // "viewDistrib*"
QT_MOC_LITERAL(10, 104, 16) // "toBeDisconnected"

    },
    "DistribWidget\0ComboChanged\0\0allowLoading\0"
    "LoadData\0fileToOpen\0update\0ShowPlot\0"
    "restoreParent\0viewDistrib*\0toBeDisconnected"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_DistribWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   44,    2, 0x0a /* Public */,
       1,    0,   47,    2, 0x2a /* Public | MethodCloned */,
       4,    1,   48,    2, 0x0a /* Public */,
       6,    0,   51,    2, 0x0a /* Public */,
       7,    0,   52,    2, 0x0a /* Public */,
       8,    1,   53,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, QMetaType::Bool,    3,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    5,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 9,   10,

       0        // eod
};

void DistribWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        DistribWidget *_t = static_cast<DistribWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->ComboChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: _t->ComboChanged(); break;
        case 2: _t->LoadData((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 3: _t->update(); break;
        case 4: _t->ShowPlot(); break;
        case 5: _t->restoreParent((*reinterpret_cast< viewDistrib*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject DistribWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_DistribWidget.data,
      qt_meta_data_DistribWidget,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *DistribWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *DistribWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_DistribWidget.stringdata0))
        return static_cast<void*>(const_cast< DistribWidget*>(this));
    return QWidget::qt_metacast(_clname);
}

int DistribWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
