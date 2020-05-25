#include "Position.h"

Position::Position():Position(0,0,0)
{
}

Position::Position(double xx, double yy, double zz):x{xx},y{yy},z{zz}
{
	r2 = xx * xx + yy * yy + zz * zz;
	this->reloadSpherical();
}

Position Position::Center()
{
	return Position(0.0,0.0,0.0);
}

Position Position::generateFromSphericalCoordinates(double r, double theta, double phi)
{
	Position p;
	p.r = r;
	p.theta = theta;
	p.phi = phi;
	p.r2 = r * r;
	p.reloadCartesian();	
	return p;
}

void Position::reloadSpherical()
{
	r2 = x * x + y * y + z * z;
	r = sqrt(r2);
	if (z == 0)
		theta = 0;
	else
		theta = atan(sqrt(x*x + y * y) / z);
	if (x == 0 && y <= 0)
		phi = 3 * PI / 2;
	else if (x == 0 && y > 0)
		phi = PI / 2;
	else
		phi = atan(y / x);

}

void Position::reloadCartesian()
{
	if (theta < 0.0 || theta >= PI || phi < 0.0 || phi >= 2 * PI)
		throw new std::exception("POSITION: Error in definnig spherical coordinates");
	r2 = r * r;
	z = r * cos(theta);
	x = r * sin(theta) * cos(phi);
	y = r * sin(theta) * sin(phi);
}

Position operator+(const Position & l, const Position & r)
{
	return Position(l.x+r.x, l.y+r.y,l.z+r.z);
}

Position operator-(const Position & l, const Position & r)
{
	return Position(l.x-r.x,l.y-r.y,l.z-r.z);
}

Position operator*(double l, const Position& r)
{
	return Position(l * r.x, l * r.y, l * r.z);
}

bool operator==(const Position & l, const Position & r)
{
	return(l.x == r.x && l.y == r.y && l.z == r.z);
}

