#include "Delone.h"
#include "S6.h"
#include "LatticeConverter.h"
#include "MaximaTools.h"
#include "ReadCellData.h"
#include "Reducer.h"
#include "ToString.h"
#include "S6M_SellingReduce.h"


void LatticeConverter::TextOutput(const std::string& label, const std::string& lattice, const LRL_Cell& cell) const {
    std::cout << label << std::endl;
    std::cout << "lattice " << lattice << std::endl;
    std::cout << "LRL_Cell  " << ToString(cell) << std::endl;
    std::cout << "G6 " << ToString(cell.Cell2V6()) << std::endl;
    std::cout << "D7  " << ToString(D7(cell.Cell2V6())) << std::endl;
    std::cout << "Scalars " << ToString(S6(cell)) << std::endl;
}

void LatticeConverter::MaximaOutput(const std::string& label, const std::string& lattice, const LRL_Cell& cell) const {
    std::cout << label << std::endl;
    std::cout << "lattice " << lattice << std::endl;
    std::cout << "LRL_Cell  " << MaximaTools::MaximaFromString(ToString(cell)) << std::endl;
    std::cout << "G6 " << MaximaTools::MaximaFromString(ToString(cell.Cell2V6())) << std::endl;
    std::cout << "D7  " << MaximaTools::MaximaFromString(ToString(D7(cell.Cell2V6()))) << std::endl;
    std::cout << "Scalars " << MaximaTools::MaximaFromString(ToString(S6(cell))) << std::endl;
}

LatticeConverter::LatticeConverter(const eOutputType type)
        : m_OutputType(type)
    {

    }

    void LatticeConverter::SetOutputMaxima(void) { m_OutputType = emaxima; }
    void LatticeConverter::SetOutputText(void) { m_OutputType = etext; }

    void LatticeConverter::Output(const std::string& label, const std::string& lattice, const LRL_Cell& cell) const {
        (m_OutputType == etext) ? TextOutput(label, lattice, cell) : MaximaOutput(label, lattice, cell);
        //std::cout << label << std::endl;
        //std::cout << "lattice " << lattice << std::endl;
        //std::cout << "LRL_Cell  " << ToString(cell) << std::endl;
        //std::cout << "G6 " << ToString(cell.Cell2V6()) << std::endl;
        //std::cout << "D7  " << ToString(D7(cell.Cell2V6())) << std::endl;
        //std::cout << "Scalars " << ToString(S6(cell)) << std::endl;
    }

    LRL_Cell LatticeConverter::NiggliReduceCell(const std::string& lattice, const LRL_Cell& cell) {
        const G6 g6 = cell.Cell2V6();
        const Mat66 mLattice = LRL_Cell(cell).LatSymMat66(lattice);
        Mat66 m66;
        G6 g6Vec;
        G6 redVec;
        int reduced;
        g6Vec=mLattice*g6;
        CS6M_G6Reduce(g6Vec,redVec,reduced);
        /* const bool b = Reducer::Reduce(mLattice*g6, m66, redVec, 0.00000001); */
        if (reduced) {
            return LRL_Cell(redVec);
        }
        else {
            return LRL_Cell();
        }
    }

    void LatticeConverter::NiggliReducedOutput(const std::string& label, const std::string& lattice, const LRL_Cell& cell) {
        const LRL_Cell reducedLRL_Cell = NiggliReduceCell(lattice, cell);
        Output(label, "P", reducedLRL_Cell);
    }

    LRL_Cell LatticeConverter::DeloneReduceCell(const std::string& lattice, const LRL_Cell& cell) {
        const G6 g6 = cell.Cell2V6();
        const Mat66 mLattice = LRL_Cell(cell).LatSymMat66(lattice);
        Mat66 m66;
        G6 g6Vec;
        G6 redVec;
        D7 d7Vec;
        D7 d7redVec;
        int reduced;
        g6Vec=mLattice*g6;
        CS6M_G6toD7(g6Vec,d7Vec);
        CS6M_D7Reduce(d7Vec,d7redVec,reduced);
        /* const bool b = Delone::Reduce(mLattice*g6, m66, redVec, 0.00000001); */
        if (reduced) {
            CS6M_D7toG6(d7redVec,redVec);
            return LRL_Cell(redVec);
        }
        else {
            return LRL_Cell();
        }
    }

    void LatticeConverter::DeloneReducedOutput(const std::string& label, const std::string& lattice, const LRL_Cell& cell) {
        const LRL_Cell reducedLRL_Cell = DeloneReduceCell(lattice, cell);
        Output(label, "P", reducedLRL_Cell);
    }
