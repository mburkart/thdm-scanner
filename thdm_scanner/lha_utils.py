#!/usr/bin/env python

import math

from thdm_scanner.utility.utils import fortran_s_to_d, d_to_fortran_s


class LHAEntry(object):
    """Class representing an entry in an LHA input or output file"""

    def __init__(self, name, value, comment=""):
        self._name = name
        self._value = value
        self._comment = comment

    @property
    def name(self):
        return self._name

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val):
        self._value = val

    def __repr__(self):
        return "LHAEntry(name={n}, value={val})".format(n=self._name,
                                                        val=self._value)

    def __str__(self):
        if "," in self.name:
            entry_str = "{name1:>5}{name2:>5}\t{val} #  {com}".format(
                            name1=self.name.split(",")[0],
                            name2=self.name.split(",")[1],
                            val=self._value,
                            com=self._comment)
        else:
            entry_str = "{name:>5}\t{val} #  {com}".format(
                            name=self.name,
                            val=self._value,
                            com=self._comment)
        return entry_str

    def __eq__(self, other):
        return self._name == other.name


class LHABlock(object):
    """Class representing a block in an LHA input or output file.

    A block consists of one or more entries.
    """

    def __init__(self, name, comment=""):
        self._name = name
        self._comment = comment
        self._entries = []

    def __repr__(self):
        return "LHABlock(name={n})".format(n=self._name)

    def __str__(self):
        if self._comment == "":
            blk_str = "Block {name}".format(
                        name=self._name)
        else:
            blk_str = "Block {name} # {com}".format(
                        name=self._name,
                        com=self._comment)

        for entry in self._entries:
            blk_str += "\n{}".format(entry)
        return blk_str

    @property
    def name(self):
        return self._name

    def add_entry(self, entry):
        if entry in self._entries:
            raise ValueError("Entry {} already in block.".format(repr(entry)))
        else:
            self._entries.append(entry)

    def add_entry_from_vals(self, *vals):
        # First create entry from given values
        entry = LHAEntry(*vals)
        # Then add it to list of entries
        self.add_entry(entry)

    def get_value(self, key):
        for ent in self._entries:
            if ent.name == key:
                return ent.value
        else:
            raise KeyError("No entry with name `{}` in {}".format(
                            key, self))

    def set_value(self, key, value):
        for ent in self._entries:
            if ent.name == key:
                ent.value = value
                break
        else:
            raise KeyError("No entry with name `{}` in {}".format(
                            key, self))


class DecayEntry(object):
    """Class representing an decay entry in an LHA output file"""

    def __init__(self, prod1, prod2, br, nda):
        self._decay_products = (prod1, prod2)
        self._br = br
        self._nda = nda

    @property
    def decay_products(self):
        return self._decay_products

    @property
    def br(self):
        return self._br

    def __repr__(self):
        return "DecayEntry(decay={dec}, BR={br})".format(
                    dec=self.decay_products,
                    br=self._br)

    def __str__(self):
        return "{br:>5}\t{nda}\t{id1}\t{id2}".format(
                    br=self._br,
                    nda=self._nda,
                    id1=self.decay_products[0],
                    id2=self.decay_products[1])

    def __eq__(self, other):
        return self.decay_products[0] == other.decay_products[0] \
                and self.decay_products[1] == other.decay_products[1]

    def is_tau_decay(self):
        # Check if decay is di-tau decay based on PDG IDs
        return self.decay_products == (15, -15) \
                or self.decay_products == (-15, 15)

    def is_symmetric_decay(self):
        # Check if decay is into the same type of particle
        return abs(self.decay_products[0]) == abs(self.decay_products[1])


class DecayBlock(object):
    """Class representing a block in an LHA input or output file.

    A block consists of one or more entries.
    """

    _inline_comment = "#\t BR \t NDA \t ID1 \t ID2"

    def __init__(self, particle, width, comment=""):
        self._particle = particle
        self._width = width
        self._comment = comment
        self._br_ratios = []

    def __repr__(self):
        return "DecayBlock(particle={p})".format(p=self._particle)

    def __str__(self):
        dec_str = "DECAY\t{part}\t{width}\t# {com}\n{icom}".format(
                    part=self._particle,
                    width=self._width,
                    com=self._comment,
                    icom=self._inline_comment)
        for br_ratio in self._br_ratios:
            dec_str += "\n{}".format(br_ratio)
        return dec_str

    def add_branching_ratio(self, br_ratio):
        if br_ratio in self._br_ratios:
            raise ValueError("Entry already in decay block.")
        else:
            self._br_ratios.append(br_ratio)

    def get_branching_ratio(self, dec_prods):
        for entry in self._br_ratios:
            if entry._decay_products == dec_prods:
                return entry.br
        else:
            print("No decay with products {} in:\n{}".format(
                            dec_prods, self))
            print("Setting BR to 0.")
            return "0."


class LHAFile(object):
    """Class representing an LHA input or output file.

    A file consists of multiple blocks and possibly decay results.
    """

    def __init__(self, filename):
        self._filename = filename
        self._blocks = []
        self._comment = ""

    def __repr__(self):
        return "LHAFile(fname={})".format(self.filename)

    @property
    def filename(self):
        return self._filename

    def read_file(self):
        with open(self.filename, "r") as fi:
            _first_blk_found = False
            _last_was_blk = None
            for line in fi:
                line = line.rstrip("\n")
                cmnt = ""
                if not (line.startswith("Block")
                        or line.startswith("DECAY")
                        or line.startswith("BLOCK")) and not _first_blk_found:
                    if self._comment == "":
                        self._comment += "{}".format(line)
                    else:
                        self._comment += "\n{}".format(line)
                elif line.lstrip().startswith("#"):
                    continue
                elif (line.startswith("Block")
                        or line.startswith("BLOCK")):
                    _first_blk_found = True
                    if "#" in line:
                        _, blk_name, cmnt = line.split(None, 2)  # for python 3 use maxsplit  # noqa: E501
                        cmnt = cmnt.lstrip("# ")
                    else:
                        _, blk_name = line.split()
                    self._blocks.append(LHABlock(blk_name, cmnt))
                    _last_was_blk = True
                elif line.startswith("DECAY"):
                    _first_blk_found = True
                    _, particle, width, cmnt = line.split(None, 3)
                    cmnt = cmnt.lstrip("# ")
                    self._blocks.append(DecayBlock(particle, width, cmnt))
                    _last_was_blk = False
                else:
                    if _last_was_blk:
                        if len(line.split("#")[0].split()) == 1:
                            name = ""
                            value = line.split("#")[0].strip()
                            cmnt = line.split("#")[1].strip()
                        elif len(line.split("#")[0].split()) == 3:
                            # Special treatment for CKM matrix element
                            # in input file
                            name1, name2, value, cmnt = line.split(None, 3)
                            name = ",".join([name1, name2])
                        else:
                            name, value, cmnt = line.split(None, 2)
                            cmnt = cmnt.lstrip("# ")
                        self._blocks[-1].add_entry(LHAEntry(name, value, cmnt))
                    else:
                        br, nda, prod1, prod2 = line.split()
                        self._blocks[-1].add_branching_ratio(
                                DecayEntry(prod1, prod2, br, nda))
        return

    def write_file(self, fname=None):
        # Allow writing to a different file than the one that was read.
        with open(self.filename if fname is None else fname, "w") as fo:
            _wrote_line = False
            if len(self._comment) > 0:
                fo.write(self._comment)
                _wrote_line = True
            for block in self._blocks:
                if _wrote_line:
                    fo.write("\n")
                fo.write(str(block))
                _wrote_line = True
        return

    def _get_entry_value(self, blk_name, key):
        for block in self._blocks:
            if isinstance(block, LHABlock):
                if block.name == blk_name:
                    return block.get_value(key)
            else:
                if block._particle == blk_name:
                    br_key = tuple(key.split(","))
                    br = block.get_branching_ratio(br_key)
                    return br
        else:
            raise KeyError("No block `{}` in {}".format(blk_name, self))

    def _set_entry_value(self, blk_name, key, value,
                         choices=None):
        if choices is not None:
            if value not in choices:
                err_mssg = ("Given value {} for key {} in block {} not allowed. "  # noqa: E501
                            "Choose from {}".format(
                                value, key,
                                blk_name, choices))
                raise ValueError(err_mssg)
        for block in self._blocks:
            if isinstance(block, LHABlock):
                if block.name == blk_name:
                    block.set_value(key, value)
                    break
            else:
                raise NotImplementedError(
                        "Setting of branching ratios not supported")
        else:
            raise KeyError("No block `{}` in {}".format(blk_name, self))


class SusHiInput(LHAFile):

    def __init__(self, filename):
        super(SusHiInput, self).__init__(filename)
        # sushi inputs
        sushi_block = LHABlock("SUSHI")
        sushi_block.add_entry_from_vals("1", "2", "model: 0 = SM, 1 = MSSM, 2 = 2HDM, 3 = NMSSM")  # noqa: E501
        sushi_block.add_entry_from_vals("2", "11", "11 = h, 12 = H, 21 = A")
        sushi_block.add_entry_from_vals("3", "0", "collider: 0 = p-p, 1 = p-pbar")  # noqa: E501
        sushi_block.add_entry_from_vals("4", "13000.d0", "center-of-mass energy in GeV")  # noqa: E501
        sushi_block.add_entry_from_vals("5", "2", "order ggh: -1 = off, 0 = LO, 1 = NLO, 2 = NNLO, 3 = N3LO")  # noqa: E501
        sushi_block.add_entry_from_vals("6", "2", "order bbh: -1 = off, 0 = LO, 1 = NLO, 2 = NNLO")  # noqa: E501
        sushi_block.add_entry_from_vals("7", "1", "electroweak cont. for ggh:")
        sushi_block.add_entry_from_vals("19", "0", "0 = silent mode of SusHi, 1 = normal output")  # noqa: E501
        sushi_block.add_entry_from_vals("20", "0", "ggh@nnlo subprocesses: 0=all, 10=ind. contributions")  # noqa: E501
        # 2HDMC input block
        thdmc_block = LHABlock("2HDMC", "2HDMC arXiv:0902.0851")
        thdmc_block.add_entry_from_vals("-1", "0", "CMD line mode: 0 direct link to library, 1 command line mode")  # noqa: E501
        thdmc_block.add_entry_from_vals("1", "3", "2HDMC key, 1=lambda basis, 2=physical basis, 3=H2 basis")  # noqa: E501
        thdmc_block.add_entry_from_vals("2", "2", "2HDM version type: (1=Type I,2=Type II,3=Flipped,4=Lepton Specific)")  # noqa: E501
        thdmc_block.add_entry_from_vals("3", "10.", "tan(beta)")
        thdmc_block.add_entry_from_vals("4", "100.", "m12")
        thdmc_block.add_entry_from_vals("21", "125.38d0", "mh")
        thdmc_block.add_entry_from_vals("22", "300d0", "mH")
        thdmc_block.add_entry_from_vals("23", "400d0", "mA")
        thdmc_block.add_entry_from_vals("24", "400d0", "mC")
        thdmc_block.add_entry_from_vals("25", "0.995", "sin(beta-alpha)")
        thdmc_block.add_entry_from_vals("26", "0.0d0", "lambda_6")
        thdmc_block.add_entry_from_vals("27", "0.0d0", "lambda_7")
        thdmc_block.add_entry_from_vals("31", "125.38d0", "mh")
        thdmc_block.add_entry_from_vals("32", "200.d0", "mH")
        thdmc_block.add_entry_from_vals("33", "0.5d0", "sin(beta-alpha)")
        thdmc_block.add_entry_from_vals("34", "0.1d0", "Z4")
        thdmc_block.add_entry_from_vals("35", "0.1d0", "Z5")
        thdmc_block.add_entry_from_vals("36", "0.1d0", "Z7")
        # SM Input block
        sm_block = LHABlock("SMINPUTS", "Standard Model inputs")
        sm_block.add_entry_from_vals("1", "1.27934000e+02", "alpha_em^(-1)(MZ) SM MSbar")  # noqa: E501
        sm_block.add_entry_from_vals("2", "1.16637000e-05", "G_Fermi")
        # sm_block.add_entry_from_vals("3", "1.17200000e-01", "alpha_s(MZ) SM MSbar")  # noqa: E501
        sm_block.add_entry_from_vals("3", "1.18000000e-01", "alpha_s(MZ) SM MSbar")  # noqa: E501
        sm_block.add_entry_from_vals("4", "9.11876000e+01", "m_Z(pole)")
        sm_block.add_entry_from_vals("5", "4.18000000e+00", "m_b(m_b)")
        sm_block.add_entry_from_vals("6", "1.72500000e+02", "m_t(pole)")
        sm_block.add_entry_from_vals("8", "1.27900000e+00", "m_c(m_c)")
        # distribution blocPDF4LHC15_nnlo_mc
        dist_block = LHABlock("DISTRIB")
        dist_block.add_entry_from_vals("1", "0", "distribution : 0 = sigma_total, 1 = dsigma/dpt,")  # noqa: E501
        dist_block.add_entry_from_vals("2", "0", "pt-cut: 0 = no, 1 = pt > ptmin, 2 = pt < ptmax,")  # noqa: E501
        dist_block.add_entry_from_vals("21", "30.d0", "minimal pt-value ptmin in GeV")  # noqa: E501
        dist_block.add_entry_from_vals("22", "100.d0", "maximal pt-value ptmax in GeV")  # noqa: E501
        dist_block.add_entry_from_vals("3", "0", "rapidity-cut: 0 = no, 1 = Abs[y] < ymax,")  # noqa: E501
        dist_block.add_entry_from_vals("31", "0.5d0", "minimal rapidity ymin")
        dist_block.add_entry_from_vals("32", "1.5d0", "maximal rapidity ymax")
        dist_block.add_entry_from_vals("4", "0", "0 = rapidity, 1 = pseudorapidity")  # noqa: E501
        # Scales input block
        scale_block = LHABlock("SCALES")
        scale_block.add_entry_from_vals("1", "0.5", "renormalization scale muR/mh")  # noqa: E501
        scale_block.add_entry_from_vals("2", "0.5", "factorization scale muF/mh")  # noqa: E501
        scale_block.add_entry_from_vals("11", "1.0", "renormalization scale muR/mh for bbh")  # noqa: E501
        scale_block.add_entry_from_vals("12", "0.25", "factorization scale muF/mh for bbh")  # noqa: E501
        scale_block.add_entry_from_vals("3", "0", "1 = Use (muR,muF)/Sqrt(mh^2+pt^2) for dsigma/dpt and d^2sigma/dy/dpt")  # noqa: E501
        # Bottom renormalization input
        renorm_block = LHABlock("RENORMBOT", "Renormalization of the bottom sector")  # noqa: E501
        renorm_block.add_entry_from_vals("1", "0", "m_b used for bottom Yukawa:  0 = OS, 1 = MSbar(m_b), 2 = MSbar(muR)")  # noqa: E501
        renorm_block.add_entry_from_vals("4", "4.75d0", "Fixed value of m_b^OS")  # noqa: E501
        # PDF input block
        pdf_block = LHABlock("PDFSPEC")
        pdf_block.add_entry_from_vals("1", "MMHT2014lo68cl.LHgrid", "name of pdf (lo)")  # noqa: E501
        pdf_block.add_entry_from_vals("2", "PDF4LHC15_nlo_mc_pdfas.LHgrid", "name of pdf (nlo)")  # noqa: E501
        pdf_block.add_entry_from_vals("3", "PDF4LHC15_nnlo_mc_pdfas.LHgrid", "name of pdf (nnlo)")  # noqa: E501
        pdf_block.add_entry_from_vals("4", "PDF4LHC15_nnlo_mc_pdfas.LHgrid", "name of pdf (n3lo)")  # noqa: E501
        # pdf_block.add_entry_from_vals("10", "0", "set number - if different for LO, NLO, NNLO, N3LO use entries 11, 12, 13")  # noqa: E501
        pdf_block.add_entry_from_vals("11", "0", "set number - if different for LO, NLO, NNLO, N3LO use entries 11, 12, 13")  # noqa: E501
        pdf_block.add_entry_from_vals("12", "0", "set number - if different for LO, NLO, NNLO, N3LO use entries 11, 12, 13")  # noqa: E501
        pdf_block.add_entry_from_vals("13", "0", "set number - if different for LO, NLO, NNLO, N3LO use entries 11, 12, 13")  # noqa: E501
        # vegas input block
        vegas_block = LHABlock("VEGAS")
        vegas_block.add_entry_from_vals("1", "10000", "Number of points")
        vegas_block.add_entry_from_vals("2", "5", "Number of iterations")
        vegas_block.add_entry_from_vals("3", "10", "Output format of VEGAS integration")  # noqa: E501
        vegas_block.add_entry_from_vals("4", "2000", "Number of points")
        vegas_block.add_entry_from_vals("5", "5", "Number of iterations")
        vegas_block.add_entry_from_vals("14", "5000", "Number of points in second run")  # noqa: E501
        vegas_block.add_entry_from_vals("15", "2", "Number of iterations in second run")  # noqa: E501
        vegas_block.add_entry_from_vals("6", "0", "Output format of VEGAS integration")  # noqa: E501
        vegas_block.add_entry_from_vals("7", "2000", "Number of points")
        vegas_block.add_entry_from_vals("8", "5", "Number of iterations")
        vegas_block.add_entry_from_vals("17", "5000", "Number of points in second run")  # noqa: E501
        vegas_block.add_entry_from_vals("18", "2", "Number of iterations in second run")  # noqa: E501
        vegas_block.add_entry_from_vals("9", "0", "Output format of VEGAS integration")  # noqa: E501
        # factors block
        fact_block = LHABlock("FACTORS")
        fact_block.add_entry_from_vals("1", "0.d0", "factor for yukawa-couplings: c")  # noqa: E501
        fact_block.add_entry_from_vals("2", "1.d0", "t")
        fact_block.add_entry_from_vals("3", "1.d0", "b")
        # block addition
        self._blocks.append(sushi_block)
        self._blocks.append(thdmc_block)
        self._blocks.append(sm_block)
        self._blocks.append(dist_block)
        self._blocks.append(scale_block)
        self._blocks.append(renorm_block)
        self._blocks.append(pdf_block)
        self._blocks.append(vegas_block)
        self._blocks.append(fact_block)

    def __repr__(self):
        return "SusHiInput(fname={})".format(self.filename)

    @property
    def higgs_boson(self):
        return int(self._get_entry_value("SUSHI", "2"))

    @higgs_boson.setter
    def higgs_boson(self, higgs_code):
        self._set_entry_value("SUSHI", "2", str(higgs_code),
                              choices={"11", "12", "21"})

    @property
    def thdm_basis(self):
        return int(self._get_entry_value("2HDMC", "1"))

    @thdm_basis.setter
    def thdm_basis(self, basis):
        self._set_entry_value("2HDMC", "1", str(basis),
                              choices={"1", "2", "3"})

    @property
    def thdm_type(self):
        return int(self._get_entry_value("2HDMC", "2"))

    @thdm_type.setter
    def thdm_type(self, thdm_type):
        self._set_entry_value("2HDMC", "2", str(thdm_type),
                              choices={"1", "2"})

    @property
    def tanb(self):
        return float(self._get_entry_value("2HDMC", "3"))

    @tanb.setter
    def tanb(self, tanb):
        self._set_entry_value("2HDMC", "3", str(tanb))

    @property
    def mh(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "31"))

    @mh.setter
    def mh(self, mh):
        self._set_entry_value("2HDMC", "31", d_to_fortran_s(mh))
        self._set_entry_value("2HDMC", "21", d_to_fortran_s(mh))

    @property
    def mH(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "32"))

    @mH.setter
    def mH(self, mH):
        self._set_entry_value("2HDMC", "32", d_to_fortran_s(mH))
        self._set_entry_value("2HDMC", "22", d_to_fortran_s(mH))

    @property
    def sin_betal(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "33"))

    @sin_betal.setter
    def sin_betal(self, sin):
        self._set_entry_value("2HDMC", "33", d_to_fortran_s(sin))
        self._set_entry_value("2HDMC", "25", d_to_fortran_s(sin))

    @property
    def Z4(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "34"))

    @Z4.setter
    def Z4(self, z4):
        self._set_entry_value("2HDMC", "34", d_to_fortran_s(z4))

    @property
    def Z5(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "35"))

    @Z5.setter
    def Z5(self, z5):
        self._set_entry_value("2HDMC", "35", d_to_fortran_s(z5))

    @property
    def Z7(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "36"))

    @Z7.setter
    def Z7(self, z7):
        self._set_entry_value("2HDMC", "36", d_to_fortran_s(z7))

    @property
    def m12(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "4"))

    @m12.setter
    def m12(self, m12):
        self._set_entry_value("2HDMC", "4", d_to_fortran_s(m12))

    @property
    def mA(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "23"))

    @mA.setter
    def mA(self, mA):
        self._set_entry_value("2HDMC", "23", d_to_fortran_s(mA))

    @property
    def mHp(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "24"))

    @mHp.setter
    def mHp(self, mHp):
        self._set_entry_value("2HDMC", "24", d_to_fortran_s(mHp))

    @property
    def lambda6(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "26"))

    @lambda6.setter
    def lambda6(self, lambda6):
        self._set_entry_value("2HDMC", "26", d_to_fortran_s(lambda6))

    @property
    def lambda7(self):
        return fortran_s_to_d(self._get_entry_value("2HDMC", "27"))

    @lambda7.setter
    def lambda7(self, lambda7):
        self._set_entry_value("2HDMC", "27", d_to_fortran_s(lambda7))
    # TODO: Check these new quantities
    # what does entry 50 do?

    @property
    def pdf_set(self):
        return self._get_entry_value("PDFSPEC", "13")

    @pdf_set.setter
    def pdf_set(self, pdf_set):
        self._set_entry_value("PDFSPEC", "12", str(pdf_set),
                              choices=map(str, range(103)))
        self._set_entry_value("PDFSPEC", "13", str(pdf_set),
                              choices=map(str, range(103)))

    @property
    def alpha_s(self):
        return self._get_entry_value("SMINPUTS", "3")

    @alpha_s.setter
    def alpha_s(self, alpha_s):
        self._set_entry_value("SMINPUTS", "3", alpha_s)


class THDMCOutput(LHAFile):

    def __init__(self, filename):
        super(THDMCOutput, self).__init__(filename)
        self.read_file()

    @property
    def valid_model(self):
        val_flags = []
        for key in range(1, 5):
            val_flags.append(int(self._get_entry_value("THDM", str(key))))
        return sum(val_flags) == 4

    @property
    def mh(self):
        return float(self._get_entry_value("MASS", "25"))

    @property
    def mH(self):
        return float(self._get_entry_value("MASS", "35"))

    @property
    def mA(self):
        return float(self._get_entry_value("MASS", "36"))

    @property
    def mHp(self):
        return float(self._get_entry_value("MASS", "37"))

    @property
    def br_htautau(self):
        return float(self._get_entry_value("25", "15,-15"))

    @property
    def br_Htautau(self):
        return float(self._get_entry_value("35", "15,-15"))

    @property
    def br_Atautau(self):
        return float(self._get_entry_value("36", "15,-15"))


class SusHiOutput(LHAFile):

    def __init__(self, filename):
        super(SusHiOutput, self).__init__(filename)
        self.read_file()

    @property
    def xs_ggPhi(self):
        return float(self._get_entry_value("SUSHIggh", "1"))

    @property
    def xs_ggPhi_scale_up(self):
        return float(self._get_entry_value("SUSHIggh", "103"))

    @property
    def xs_ggPhi_scale_down(self):
        return float(self._get_entry_value("SUSHIggh", "102"))

    @property
    def xs_bbPhi(self):
        return float(self._get_entry_value("SUSHIbbh", "1"))

    @property
    def mPhi(self):
        return float(self._get_entry_value("MASSOUT", "1"))
