/*
 * File: gzIO.cpp
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include "gzIO.h"
#include <iostream>

bool readline(gzFile f, string& l, int len)
{
    // std::vector< char > v(256);
    string      v(len, ' ');
    std::size_t pos = 0;

    if (gzgets(f, &v[pos], v.size()) == 0)
    {
        // end-of-file or error
        int         err;
        const char* msg = gzerror(f, &err);
        if (err != Z_OK)
        {
            cout << "read gz file error, error_code: " << err << " error_msg: " << msg << endl;
            return false;
        }
    }

    pos = v.find('\n');
    if (pos == string::npos)
        pos = 0;

    v.resize(pos);
    l = v;

    v.clear();
    if (l.size() != 0)
        return true;
    return false;
}
