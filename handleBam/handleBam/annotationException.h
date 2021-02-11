/*
 * File: annotationException.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <exception>
#include <string>

class AnnotationException : public std::logic_error
{
public:
    AnnotationException(const std::string& s) : std::logic_error(s) {}
};