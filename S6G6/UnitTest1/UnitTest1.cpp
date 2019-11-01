#include "pch.h"
#include "CppUnitTest.h"
#include "LatticeConverter.h"
#include "CS6Dist.h"
#include "S6.h"
#include <fstream>
#include <iostream>
#include <S6Dist.h>
#include <windows.h>



using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{
	TEST_CLASS(UnitTest1)
	{
	public:

		//bounds for which cells to test from PDB cell database
		int min = 1;
		int max = 3;

		void test_distances(const std::vector<S6> cells, const std::vector<double> expectedDistances) {
			int c = 0;
			for (int i = 0; i < cells.size(); i++) {
				for (int j = i + 1; j < cells.size(); j++) {
					double dist = CS6Dist(cells[i].data(), cells[j].data());
					Assert::AreEqual(round(expectedDistances[c] * 100) / 100, round(dist * 100) / 100);
					c++;
				}
			}
		}

		void testS6CS6(const std::vector<S6> cells) {
			S6Dist s6dist(0);
			int n_discrepancies = 0;
			std::vector<double> errors;
			std::vector<S6> cells1;
			std::vector<S6> cells2;

			for (int i = 0; i < cells.size(); i++) {
				for (int j = i + 1; j < cells.size(); j++) {
					double cs6distance = CS6Dist(cells[i].data(), cells[j].data());
					double s6distance = s6dist.DistanceBetween(cells[i], cells[j]);
					if (round(cs6distance * 100) / 100 != round(s6distance * 100) / 100) {
						n_discrepancies++;
						cells1.push_back(cells[i]);
						cells2.push_back(cells[j]);
						errors.push_back(cs6distance - s6distance);
					}
					//Assert::AreEqual(round(cs6distance * 100) / 100, round(s6distance * 100) / 100);
					
				}
			}

			Assert::AreEqual(0, n_discrepancies);
		}
		TEST_METHOD(centeringSymbols)
		{
			LatticeConverter converter;
			const std::vector<std::string> centeringSymbols = { "p", "f", "i" };
			const std::vector<std::string> cells = { "10 10 10 90 90 90","10 10 10 90 90 90" ,"10 10 10 90 90 90" };
			const std::vector<double> expectedDistances = { 132.288, 136.931, 35.3553 };
			std::vector<S6> reducedCells;
			for (int i = 0; i < centeringSymbols.size(); i++) {
				LRL_Cell cell = converter.SellingReduceCell(centeringSymbols[i], cells[i]);
				S6 cellS6(cell);
				reducedCells.push_back(cellS6);
			}
			test_distances(reducedCells, expectedDistances);
		}
		TEST_METHOD(CS6vsS6)
		{
			LatticeConverter converter;
			const std::vector<std::string> cells = {
				"p 127.910 87.420 130.270 90 119.49 90",
				"p 75.873 46.814 83.678 90 90.93 90",
				"p 75.340 75.340 120.610 90 90 120",
				"p 93.910 93.910 131.170 90 90 120",
				"p 94.145 94.145 131.616 90 90 120",
				"p 92.520 92.520 128.946 90 90 120",
				"p 92.787 92.787 130.655 90 90 120",
				"p 92.718 92.718 129.965 90 90 120",
				"p 180.868 180.868 177.111 90 90 90",
				"p 48.150 81.540 78.500 90 104.75 90",
				"c 173.510 62.020 125.670 90 99.92 90",
				"p 51.528 104.795 164.389 90 90 90",
				"p 55.510 83.270 111.040 90 90 90",
				"i 687.9 687.9 1933 90 90 90",
				"p 334.67 597.13 336.5 90 120.29 90"
			};
			std::vector<S6> reducedCells;
			for (std::string buffer : cells) {
				int i = buffer.find_first_of(" ");
				std::string centeringSymbol = buffer.substr(0, i);
				std::string cell = buffer.substr(i);
				LRL_Cell lrl = converter.SellingReduceCell(centeringSymbol, cell);
				S6 cellS6(lrl);
				reducedCells.push_back(cellS6);
			}
			testS6CS6(reducedCells);
		}

		TEST_METHOD(readCSV) {
			std::string filename = "..\\PDBcelldatabase.csv";
			std::ifstream in(filename);
			LatticeConverter converter;

			std::vector <std::vector <std::string> > cells;
			int row = 0;
			while (in)
			{
				std::string s;
				if (!getline(in, s)) break;

				std::istringstream ss(s);
				std::vector <std::string> record;

				int col = 0;
				while (ss)
				{
					std::string token;
					if (!getline(ss, token, ',')) break;
					if (col > 0 && col < 8 && row > 0) {
						token = token.substr(1, token.size()-2);
						record.push_back(token.c_str());
					}
					col++;
				}

				if(row >= min && row < max) cells.push_back(record);
				row++;
			}
			//Assert::IsTrue(in.eof());

			std::vector<S6> reducedCells;
			for (std::vector<std::string> cell : cells) {
				std::string centeringSymbol(cell[6].substr(0,1));
				centeringSymbol[0] += 32; //convert to lowercase
				std::ostringstream cell_stream;
				cell_stream << " " << cell[3] << " " << cell[4] << " " << cell[5] << " " << cell[0] << " " << cell[1] << " " << cell[2];
				LRL_Cell lrl = converter.SellingReduceCell(centeringSymbol, cell_stream.str());
				S6 cellS6(lrl);
				reducedCells.push_back(cellS6);
			}
			testS6CS6(reducedCells);
		}

	};
}
