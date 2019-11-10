// auto save input history to file and repeat input when needed

#pragma once
#include "global.h"
#include "time.h"

// read a variable from stdin and print to std::out
#define SLS_INP(var) sls_inp_imp(#var, var);

namespace slisc {

// implementation for SLS_CIN
template <class T>
void sls_inp_imp(const Char * varname, T &var)
{
    cin >> var;
    cout << varname << " = " << var << ";    ";
}

// jump to the next line
inline void cin_line()
{
	cin.ignore(1000, '\n');
    cout << endl;
}

class Input
{
private:
    // file IO status
    // 0: no file IO
    // 1: reading from file
    // 2: writing to file
    enum class Stat: Char { NO_IO, WRITE, READ };
    Stat m_status;
    Str m_fname;
    ofstream m_fout;
    ifstream m_fin;

public:

    // default: no file IO
    Input() : m_status(Stat::NO_IO)
    {}

    // create new output file
    // overwrite if exists
    void newfile(Str_I fname)
    {
        if (m_status == Stat::READ) {
            m_fin.close();
        }
        else if (m_status == Stat::WRITE) {
            m_fout.close();
        }
        m_fname = fname;
        m_fout.open(fname);
        m_status = Stat::WRITE;
    }

    // open existing file for read
    // create new if not exist
    // will continue writing if more input required
    void openfile(Str_I fname)
    {
        if (m_status == Stat::READ) {
            m_fin.close();
        }
        else if (m_status == Stat::WRITE) {
            m_fout.close();
        }
        m_fname = fname;
        m_fin.open(fname);
        m_status = Stat::READ;
    }


    void Bool_cin(Bool_O out, Str_I prompt) {
        using namespace std;
        string str;
        Int i;
        for (i = 0; i < 10000; ++i) {
            cout << prompt << " (y/n) ";
            getline(cin, str);
            if (str.size() == 1) {
                if (str[0] == 'y') {
                    out = true;
                    return;
                }
                else if (str[0] == 'n') {
                    out = false;
                    return;
                }
            }
            // failed
            SLS_WARN("illegal input, try again!");
            pause(1);
        }
    }

    // return 0 if successful
    // return -1 reached empty line, m_fin position unchanged
    Int Bool_fin(Bool_O out, Str_I prompt) {
        using namespace std;
        string str;
        getline(m_fin, str);
        if (str.size() == 1) {
            if (str[0] == 'y') {
                cout << prompt << " (y/n) " << 'y' << endl;
                out = true;
                return 0;
            }
            else if (str[0] == 'n') {
                cout << prompt << " (y/n) " << 'n' << endl;
                out = false;
                return 0;
            }
        }
        // failed
        out = false;
        return -1;
    }

    void Bool_fout(Bool_O out, Str_I prompt) {
        Bool_cin(out, prompt);
        if (out)
            m_fout << 'y' << endl;
        else
            m_fout << 'n' << endl;
    }

    // input a bool
    // " (y/n) " will be appended to prompt
    Bool iBool(Str_I prompt) {
        slisc::Bool out;
        if (m_status == Stat::NO_IO) {
            // no IO
            Bool_cin(out, prompt);
        }
        else if (m_status == Stat::WRITE) {
            // write
            Bool_fout(out, prompt);
        }
        else if (m_status == Stat::READ) {
            // read
            if (Bool_fin(out, prompt) != 0) {
                m_status = Stat::WRITE;
                m_fout.open(m_fname, std::ios::app);
                m_fin.close();
                Bool_fout(out, prompt);
                return out;
            }
        }
        return out;
    }

    template <class T1, class T2>
    void num2_cin(T1 &out1, T2 &out2, Str_I prompt) {
        using namespace std;
        string str;
        stringstream ss;
        Int i;
        for (i = 0; i < 10000; ++i) {
            cout << prompt << " ";
            getline(cin, str);
            if (str.size() >= 3) {
                ss.str(str);
                ss >> out1;
                ss >> out2;
                return;
            }
            // failed
            SLS_WARN("illegal input, try again!");
            pause(1);
        }
    }

    // return 0 if successful
    // return -1 if fail
    template <class T1, class T2>
    Int num2_fin(T1 &out1, T2 &out2, Str_I prompt) {
        using namespace std;
        string str;
        stringstream ss;
        getline(m_fin, str);
        if (str.size() >= 3) {
            ss.str(str);
            ss >> out1;
            ss >> out2;
            cout << prompt << " " << out1 << " " << out2 << endl;
            return 0;
        }
        // failed
        return -1;
    }

    // return 0 if successful
    // return -1 if fail
    template <class T1, class T2>
    void num2_fout(T1 &out1, T2 &out2, Str_I prompt) {
        num2_cin(out1, out2, prompt);
        m_fout << out1 << " " << out2 << endl;
    }

    template <class T1, class T2>
    void num2(T1 &out1, T2 &out2, Str_I prompt) {
        if (m_status == Stat::NO_IO) {
            // no IO
            num2_cin(out1, out2, prompt);
        }
        else if (m_status == Stat::WRITE) {
            // write
            num2_fout(out1, out2, prompt);
        }
        else if (m_status == Stat::READ) {
            // read
            if (num2_fin(out1, out2, prompt) != 0) {
                m_status = Stat::WRITE;
                m_fout.open(m_fname, std::ios::app);
                m_fin.close();
                num2_fout(out1, out2, prompt);
            }
        }
        return;
    }
};

} // namespace slisc
