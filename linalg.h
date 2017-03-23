//
// Created by lurker on 3/22/17.
//

#ifndef LEVELSET_LINALG_H
#define LEVELSET_LINALG_H

#include "utils.h"

namespace bbfmm {
    class Vector {
    public:
        /*
         * data members
         */
        index_t _row;
        bool_t _ownership;
        scalar_t *_data;

        Vector(int row = 0) : _row(row), _ownership(true) {
            if (_row > 0) {
                _data = new scalar_t[_row];
                assert(_data != nullptr);
                memset(_data, 0, _row * sizeof(scalar_t));
            }
            else { _data = nullptr; }
        }

        Vector(int row, bool_t ownership, scalar_t *data) : _row(row), _ownership(ownership) {
            if (_ownership) {
                if (_row > 0) {
                    _data = new scalar_t[_row];
                    assert(_data != nullptr);
                    memcpy(_data, data, _row * sizeof(scalar_t));
                }
                else {
                    _data = nullptr;
                }
            } else {
                /*
                 * has to assure data holds exact "row" elements.
                 */
                _data = data;
            }
        }


        Vector(const Vector &v) : _row(v._row), _ownership(v._ownership) {
            if (_ownership) {
                if (_row > 0) {
                    _data = new scalar_t[_row];
                    assert(_data != nullptr);
                    memcpy(_data, v._data, _row * sizeof(scalar_t));
                }
                else {
                    _data = nullptr;
                }
            } else {
                _data = v._data;
            }
        }

        ~Vector() {
            if (_ownership) {
                if (_row > 0) {
                    delete[] _data;
                    _data = nullptr;
                }
            }
        }

        Vector &operator=(const Vector &v) {
            if (_ownership) {
                if (_row > 0) {
                    delete[] _data;
                    _data = nullptr;
                }
            }

            _row = v._row;
            _ownership = v._ownership;

            if (_ownership) {
                if (_row > 0) {
                    _data = new scalar_t[_row];
                    assert(_data != nullptr);
                    memcpy(_data, v._data, _row * sizeof(scalar_t));
                }
                else {
                    _data = nullptr;
                }
            } else {
                _data = v._data;
            }
            return *this;
        }

        void resize(int row) {
            assert(_ownership);
            if (row != _row) {
                if (_row > 0) {
                    delete[] _data;
                    _data = nullptr;
                }
                _row = row;
                if (_row > 0) {
                    _data = new scalar_t[_row];
                    assert(_data != nullptr);
                    memset(_data, 0, _row * sizeof(scalar_t));
                }
                else { _data = nullptr; }
            }
        }

        const scalar_t operator()(int i) const {
            assert(i >= 0 && i < _row);
            return _data[i];
        }

        scalar_t &operator()(int i) {
            assert(i >= 0 && i < _row);
            return _data[i];
        }

        scalar_t *data() {
            return _data;
        }

        int row() {
            return _row;
        }
    };


    inline std::ostream &operator<<(std::ostream &os, Vector &v) {
        os << "rows: " << v._row << std::endl;
        os.setf(std::ios_base::scientific, std::ios_base::floatfield);
        for (int i = 0; i < v._row; ++i) {
            os << " " << v(i);
        }
        os << std::endl;
        return os;
    }

    inline void setValue(Vector &v, scalar_t val) {
        std::fill_n(v.data(), v._row, val);
    }

    inline void clear(Vector &v) {
        memset(v._data, 0, v._row * sizeof(scalar_t));
    }


/*
 * column majored matrix
 */
    class Matrix {
    public:
        int _row;
        int _col;
        bool _ownership;
        scalar_t *_data;

        Matrix(int row = 0, int col = 0) : _row(row), _col(col), _ownership(true) {
            if (_row > 0 && _col > 0) {
                _data = new scalar_t[_row * _col];
                assert(_data != nullptr);
                memset(_data, 0, _row * _col * sizeof(scalar_t));
            } else { _data = nullptr; }
        }

        Matrix(int row, int col, bool_t ownership, scalar_t *data) : _row(row), _col(col), _ownership(ownership) {
            if (_ownership) {
                if (_row > 0 && _col > 0) {
                    _data = new scalar_t[_row * _col];
                    assert(_data != nullptr);
                    memcpy(_data, data, _row * _col * sizeof(scalar_t));
                } else {
                    _data = nullptr;
                }
            } else {
                /*
                 * has to assure data holds exact "row" elements.
                 */
                _data = data;
            }
        }


        Matrix(const Matrix &v) : _row(v._row), _col(v._col), _ownership(v._ownership) {
            if (_ownership) {
                if (_row > 0 && _col > 0) {
                    _data = new scalar_t[_row * _col];
                    assert(_data != nullptr);
                    memcpy(_data, v._data, _row * _col * sizeof(scalar_t));
                } else {
                    _data = nullptr;
                }
            } else {
                _data = v._data;
            }
        }

        ~Matrix() {
            if (_ownership) {
                if (_row > 0 && _col > 0) {
                    delete[] _data;
                    _data = nullptr;
                }
            }
        }

        Matrix &operator=(const Matrix &v) {
            if (_ownership) {
                if (_row > 0 && _col > 0) {
                    delete[] _data;
                    _data = nullptr;
                }
            }

            _row = v._row;
            _col = v._col;
            _ownership = v._ownership;

            if (_ownership) {
                if (_row > 0 && _col > 0) {
                    _data = new scalar_t[_row * _col];
                    assert(_data != nullptr);
                    memcpy(_data, v._data, _row * _col * sizeof(scalar_t));
                } else {
                    _data = nullptr;
                }
            } else {
                _data = v._data;
            }
            return *this;
        }


        void resize(int row, int col) {
            assert(_ownership);
            if (row != _row || col != _col) {
                if (_row > 0 && _col > 0) {
                    delete[] _data;
                    _data = nullptr;
                }
                _row = row;
                _col = col;
                if (_row > 0 && _col > 0) {
                    _data = new scalar_t[_row * _col];
                    assert(_data != nullptr);
                    memset(_data, 0, _row * _col * sizeof(scalar_t));
                } else { _data = nullptr; }
            }
        }

        const scalar_t operator()(int i, int j) const {
            assert(i >= 0 && i < _row);
            assert(j >= 0 && j <= _col);
            return _data[i + j * _row];
        }

        scalar_t &operator()(int i, int j) {
            assert(i >= 0 && i < _row);
            assert(j >= 0 && j < _col);
            return _data[i + j * _row];
        }

        scalar_t *data() { return _data; }

        scalar_t *column(int j) {
            assert(j >= 0 && j < _col);
            return &(_data[j * _row]);
        }

        int row() { return _row; }

        int col() { return _col; }

        void setColumn(int j, Vector &v) {
            assert(j >= 0 && j < _col);
            assert(v.row() == _row);
            memcpy(&(_data[j * _row]), v.data(), sizeof(scalar_t) * _row);
        }
    };


    inline std::ostream &operator<<(std::ostream &os, Matrix &v) {
        os << "rows: " << v._row << " cols: " << v._col << std::endl;
        os.setf(std::ios_base::scientific, std::ios_base::floatfield);
        for (int i = 0; i < v._row; ++i) {
            for (int j = 0; j < v._col; ++j) {
                os << " " << v(i, j);
            }
            os << std::endl;
        }
        return os;
    }

    inline void setValue(Matrix &v, scalar_t val) {
        std::fill_n(v.data(), v._row * v._col, val);
    }

    inline void clear(Matrix &v) {
        memset(v._data, 0, v._row * v._col * sizeof(scalar_t));
    }

    inline void setBlock(Matrix &A, Matrix &B, int r, int c, int rs, int cs) {
        for (int i = 0; i < rs; ++i) {
            for (int j = 0; j < cs; ++j) {
                A(i, j) = B(r + i, c + j);
            }
        }
    }
}




#endif //LEVELSET_LINALG_H

