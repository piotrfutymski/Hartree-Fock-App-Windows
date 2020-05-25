#pragma once
#include <cmath>
#include "Constatns.h"
#include <exception>

struct Position
{
	Position();

	Position(double x, double y, double z);

	Position Center();

	Position generateFromSphericalCoordinates(double r, double theta, double phi);

	void reloadSpherical();

	void reloadCartesian();

	///

	double x, y, z, r, theta, phi;
	double r2;

};

Position operator+(const Position & l, const Position & r);
Position operator-(const Position & l, const Position & r);
Position operator*(double l, const Position& r);
bool operator==(const Position & l, const Position & r);