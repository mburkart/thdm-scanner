#!/usr/bin/env python


def d_to_fortran_s(value):
    """Converts python float value to fortran real string"""
    return trim_zeros_from_dstr(str(value)) + "d0"


def fortran_s_to_d(value):
    """Converts fortran real string to python float"""
    val, exp = value.split("d")
    return float(val) * 10 ** int(exp)


def trim_zeros_from_dstr(in_str):
    int_val, dec_val = in_str.split(".")
    val_str = ".".join([int_val, dec_val.rstrip("0")])
    return val_str
