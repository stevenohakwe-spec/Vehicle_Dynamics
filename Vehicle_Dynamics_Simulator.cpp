#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <limits>
#include <SFML/Graphics.hpp>

using namespace std;

const double deltaTime = 1.00/60.00;
const double airDensity = 1.225;
constexpr double pi = 3.14159265358979323846;
//const double frictionCoefficient = 1.03;

struct SplineData
{
    double value = 0.0;
    double Mx = 0.0;
    double My = 0.0;
};

struct VehicleAccelerations
{
    double longitudinal = 0.0;  // m/s^2
    double lateral = 0.0;       // m/s^2
    double yaw = 0.0;           // rad/s^2
};

class CubicSpline
{
    private:
    
 //for (int j=0; j<)   
        vector<double> xValues;
        vector<double> yValues;
        vector<double> M;
        
        
                vector<double> naturalSecondDerivatives( const vector<double>& x, const vector<double>& y) const
                {
                    int n = static_cast<int>(x.size());
                
                    if (n != static_cast<int>(y.size()) || n < 2)
                        throw invalid_argument("Invalid spline data.");
                
                    vector<double> M(n, 0.0);
                
                    if (n == 2)
                        return M;
                
                    vector<double> lower(n, 0.0);
                    vector<double> diagonal(n, 0.0);
                    vector<double> upper(n, 0.0);
                    vector<double> rhs(n, 0.0);
                
                    diagonal[0] = 1.0;
                    diagonal[n - 1] = 1.0;
                
                    for (int i = 1; i < n - 1; ++i)
                    {
                        double h0 = x[i] - x[i - 1];
                        double h1 = x[i + 1] - x[i];
                
                        if (h0 <= 0.0 || h1 <= 0.0)
                            throw invalid_argument("Spline coordinates must increase.");
                
                        lower[i] = h0;
                        diagonal[i] = 2.0 * (h0 + h1);
                        upper[i] = h1;
                
                        rhs[i] = 6.0 * ((y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0);
                    }
                
                    // Forward elimination
                    for (int i = 1; i < n; ++i)
                    {
                        double factor = lower[i] / diagonal[i - 1];
                
                        diagonal[i] -= factor * upper[i - 1];
                        rhs[i] -= factor * rhs[i - 1];
                    }
                
                    // Back substitution
                    for (int i = n - 2; i >= 0; --i)
                    {
                        M[i] = (rhs[i] - upper[i] * M[i + 1]) / diagonal[i];
                    }
                    return M;
                }

            
                
                void calculateSecondDerivatives(vector<vector<SplineData>>& dataArray)
                {
                    if (dataArray.empty() || dataArray[0].size() < 2)
                        throw invalid_argument("Invalid spline table.");
                
                    int rows = static_cast<int>(dataArray.size());
                    int cols = static_cast<int>(dataArray[0].size());
                    
                    for (const auto& row : dataArray)
                    {
                        if (row.size() != static_cast<size_t>(cols))
                        {
                            throw invalid_argument(
                                "Spline table contains an incomplete row."
                            );
                        }
                    }
                
                    if (rows < 3)
                        throw invalid_argument("Insufficient slip-angle data.");
                
                    // Calculate vertical second derivatives (My).
                    for (int col = 1; col < cols; ++col)
                    {
                        vector<double> x;
                        vector<double> y;
                
                        for (int row = 1; row < rows; ++row)
                        {
                            x.push_back(dataArray[row][0].value);
                            y.push_back(dataArray[row][col].value);
                        }
                
                        vector<double> derivatives =
                            naturalSecondDerivatives(x, y);
                
                        for (int row = 1; row < rows; ++row)
                            dataArray[row][col].My = derivatives[row - 1];
                    }
                
                    // Calculate horizontal second derivatives (Mx).
                    for (int row = 1; row < rows; ++row)
                    {
                        vector<double> x;
                        vector<double> y;
                
                        for (int col = 1; col < cols; ++col)
                        {
                            x.push_back(dataArray[0][col].value);
                            y.push_back(dataArray[row][col].value);
                        }
                
                        vector<double> derivatives =
                            naturalSecondDerivatives(x, y);
                
                        for (int col = 1; col < cols; ++col)
                            dataArray[row][col].Mx = derivatives[col - 1];
                    }
                }


    public:
        
        CubicSpline(vector<vector<SplineData>>& dataArray) 
        {
             calculateSecondDerivatives(dataArray);
        }
        
        
    double interpolate(double entry, double headerEntry, vector<vector<SplineData>>& dataArray) const
    {
        int dataAmount = static_cast<int>(dataArray.size());
        int dataCount = static_cast<int>(dataArray[0].size());
    
        // Check whether sufficient tire data exists.
        if (dataAmount < 3 || dataCount < 2)
            throw std::invalid_argument("Insufficient tire data.");
    
        // Check whether the requested values are within the dataset.
        if (entry < dataArray[1][0].value || entry > dataArray[dataAmount - 1][0].value || headerEntry < dataArray[0][1].value || headerEntry > dataArray[0][dataCount - 1].value)
            throw std::out_of_range( "Interpolation request outside tire data range.");
    
        // STEP 1: Find the two surrounding slip-angle rows.
        int rowIndex = -1;
    
        for (int j = 1; j < dataAmount - 1; ++j)
        {
            if (entry >= dataArray[j][0].value && entry <= dataArray[j + 1][0].value)
            {
                rowIndex = j;
                break;
            }
        }
    
        if (rowIndex == -1)
            throw std::out_of_range("Could not find slip-angle interval.");
    
        // STEP 2: Perform vertical cubic interpolation.
        auto interpolateVertical = [&](int column) -> double
        {
            int j = rowIndex;
    
            double h = dataArray[j + 1][0].value - dataArray[j][0].value;
            double dx = entry - dataArray[j][0].value;
            double d = dataArray[j][column].value;
            double a = (dataArray[j + 1][column].My - dataArray[j][column].My) / (6.0 * h);
            double b = dataArray[j][column].My / 2.0;
            double c = (dataArray[j + 1][column].value - dataArray[j][column].value) / h - h * (2.0 * dataArray[j][column].My + dataArray[j + 1][column].My) / 6.0;
    
            return a * dx * dx * dx + b * dx * dx + c * dx + d;
        };
    
        // STEP 3: Interpolate vertically at EVERY load column.
        vector<double> horizontalLoads;
        vector<double> horizontalForces;
    
        for (int i = 1; i < dataCount; ++i)
        {
            horizontalLoads.push_back(dataArray[0][i].value);
            horizontalForces.push_back(interpolateVertical(i));
        }
    
        // STEP 4: Print the vertically interpolated values.
        
        if (horizontalLoads.size() == 1)
            return horizontalForces[0];
        
        vector<double> horizontalM = naturalSecondDerivatives(horizontalLoads, horizontalForces);
        
        // Handle datasets with only one load column.
    
        // STEP 5: If load exactly matches a header, return that value.
        for (int i = 1; i < dataCount; ++i)
        {
            if (headerEntry == dataArray[0][i].value)
                return horizontalForces[i - 1];
                
        }
    
    // STEP 6: Cubic interpolation across vertical load.
    for (int i = 0; i < static_cast<int>(horizontalLoads.size()) - 1; ++i)
    {
        double x0 = horizontalLoads[i];
        double x1 = horizontalLoads[i + 1];

        if (headerEntry >= x0 && headerEntry <= x1)
        {
            double h = x1 - x0;
            double dx = headerEntry - x0;

            double d = horizontalForces[i];
            double a = (horizontalM[i + 1] - horizontalM[i]) / (6.0 * h);
            double b = horizontalM[i] / 2.0;
            double c = (horizontalForces[i + 1] - horizontalForces[i]) / h - h * (2.0 * horizontalM[i] + horizontalM[i + 1]) / 6.0;

            return a * dx * dx * dx + b * dx * dx + c * dx + d;
        }
    }

    throw std::out_of_range("Could not find load interval.");

    }

};

struct DriverInput
{
    bool gasPressed = false;
    double throttle = 0.0;
    bool right = false;
    bool left = false;
    bool brakePressed = false;
    double brake = 0.0;
    double steer = 0.0;
};

struct VehicleState
{
   double positionX = 0.0;
   double positionY = 0.0;
   double heading = 0.0;
   double longitudinalVelocity = 0.0;
   double lateralVelocity = 0.0;
   double yawRate = 0.0;
   double dragForce = 0.0;
   double downForce = 0.0;
};

struct VehicleParameters
{
    double mass = 0.0;
    double wheelbase = 0.0;
    //double frontAxleDistance = 0.0;
    //double rearAxleDistance = 0.0;
    double yawInertia = 0.0;
    double downForceCoefficient = 0.0;
    double downForceFrontPercentage = 0.0;
    double dragCoefficient = 0.0;
    //double liftCoefficient = 0.0;
    double frontalArea = 0.0;
    double liftSurfaceArea = 0.0;
    double frontTrackWidth = 0.0;
    double rearTrackWidth = 0.0;
    double gravityCenterOffset[2] = {0.0};  //The first array position is center of gravity horizontal position. The second array index is center of gravity height. Positive values move the gravity center up and forward.
    double steeringLockAngle = 30; //Place holder value until the info is read from a file.
};

struct TireParameters
{
    double lateralFriction = 0.0;
    double longitudinalFriction = 0.0;
    double sideWallStiffness = 0.0;
    double sideWallHeight = 0.0;
    double treadEffectiveness = 0.0;
};

struct TireConditions
{
    double frontLeft = 0.0;
    double frontRight = 0.0;
    double rearLeft = 0.0;
    double rearRight = 0.0;
};

struct SimulationDiagnostics
{
    TireConditions verticalLoads;
    TireConditions slipAngles;
    TireConditions lateralForces;
    VehicleAccelerations accelerations;
    double lateralAcceleration = 0.0;
};

vector<vector<SplineData>> dataFile(const string& input_string)
{
    ifstream file(input_string + ".txt");

        if (!file)
            throw runtime_error("Error opening file!");
    
        string line;
    
        if (!getline(file, line))
            throw runtime_error("File is empty!");
    
        stringstream hss(line);
        string header;
    
        vector<SplineData> headerLoads;

    // Placeholder for the slip-angle column
    headerLoads.push_back(SplineData{0.0});

    // Skip "slip_angle"
    hss >> header;

    // Extract tire loads
        while (hss >> header)
        {
            size_t position = header.find('@');
    
            if (position == string::npos)
                throw runtime_error("Invalid tire load header!");
    
            string loadText = header.substr(position + 1);
            double load = stod(loadText);
    
            headerLoads.push_back(SplineData{load});
        }

    vector<vector<SplineData>> dataArray;
    dataArray.push_back(headerLoads);

    // Read numerical measurements
        while (getline(file, line))
        {
            stringstream ss(line);
            vector<SplineData> row;
            double value;
    
            while (ss >> value)
            {
                row.push_back(SplineData{value});
            }
    
            if (!row.empty())
            {
                if (row.size() != headerLoads.size())
                    throw runtime_error("Incorrect number of columns!");
    
                dataArray.push_back(row);
            }
        }

    return dataArray;
}


TireConditions SlipAngles (const VehicleParameters& currentCar,const VehicleState& car, double steeringAngle)
{
    TireConditions slipAngles;
    
    double velocityFrontLeftX = car.longitudinalVelocity - (car.yawRate * (currentCar.frontTrackWidth/2));
    double velocityFrontLeftY = car.lateralVelocity + (car.yawRate * ((currentCar.wheelbase / 2) - currentCar.gravityCenterOffset[0]));
    
    double velocityFrontRightX = car.longitudinalVelocity + (car.yawRate * (currentCar.frontTrackWidth/2));
    double velocityFrontRightY = car.lateralVelocity + (car.yawRate * ((currentCar.wheelbase / 2) - currentCar.gravityCenterOffset[0]));
    
    double velocityRearLeftX = car.longitudinalVelocity - (car.yawRate * (currentCar.rearTrackWidth/2));
    double velocityRearLeftY = car.lateralVelocity - (car.yawRate * ((currentCar.wheelbase / 2) + currentCar.gravityCenterOffset[0]));
    
    double velocityRearRightX = car.longitudinalVelocity + (car.yawRate * (currentCar.rearTrackWidth/2));
    double velocityRearRightY = car.lateralVelocity - (car.yawRate * ((currentCar.wheelbase / 2) + currentCar.gravityCenterOffset[0]));
    
    double directionFrontLeft = atan2(velocityFrontLeftY, velocityFrontLeftX) * (180/pi);
    double directionFrontRight = atan2(velocityFrontRightY, velocityFrontRightX) * (180/pi);
    double directionRearLeft = atan2(velocityRearLeftY, velocityRearLeftX) * (180/pi);
    double directionRearRight = atan2(velocityRearRightY, velocityRearRightX) * (180/pi);
    
    slipAngles.frontLeft = steeringAngle - directionFrontLeft;
    slipAngles.frontRight = steeringAngle - directionFrontRight;
    slipAngles.rearLeft = 0 - directionRearLeft;
    slipAngles.rearRight = 0 - directionRearRight;
    
    return slipAngles;
}

TireConditions calculateVerticalForces (const VehicleParameters& currentCar, VehicleState& car)
{
    TireConditions verticalLoads;
    
    double weightForceFrontAxle = (currentCar.mass * 9.81 * ((currentCar.wheelbase/2) + currentCar.gravityCenterOffset[0]))/currentCar.wheelbase;
    double weightForceRearAxle = (currentCar.mass * 9.81 * ((currentCar.wheelbase/2) - currentCar.gravityCenterOffset[0]))/currentCar.wheelbase;

    double downForceFrontAxle = (car.downForce * currentCar.downForceFrontPercentage);
    double downForceRearAxle = (car.downForce * (1- currentCar.downForceFrontPercentage));
    
    verticalLoads.frontLeft = (weightForceFrontAxle / 2) + (downForceFrontAxle / 2);
    verticalLoads.frontRight = (weightForceFrontAxle / 2) + (downForceFrontAxle / 2);
    verticalLoads.rearLeft = (weightForceRearAxle / 2) + (downForceRearAxle / 2);
    verticalLoads.rearRight = (weightForceRearAxle / 2) + (downForceRearAxle / 2);
    
    return verticalLoads;
}



TireConditions calculateLateralForces(const TireConditions& slipAngles, const TireConditions& verticalLoads, const CubicSpline& tireSpline, vector<vector<SplineData>>& tireData)
{
    TireConditions lateralForces;

    auto calculateFy = [&](double slipAngle, double verticalLoad)
    {
        double forceMagnitude = tireSpline.interpolate(std::abs(slipAngle), verticalLoad / 1000.0, tireData);

        return std::copysign(forceMagnitude, slipAngle);
    };

    lateralForces.frontLeft = calculateFy(slipAngles.frontLeft, verticalLoads.frontLeft);
    lateralForces.frontRight = calculateFy(slipAngles.frontRight, verticalLoads.frontRight);
    lateralForces.rearLeft = calculateFy(slipAngles.rearLeft, verticalLoads.rearLeft);
    lateralForces.rearRight = calculateFy(slipAngles.rearRight, verticalLoads.rearRight);

    return lateralForces;
}


TireConditions calculateMaxLongitudinalForces( const TireConditions& lateralForces, const TireConditions& verticalLoads, double longitudinalFriction, double lateralFriction)
{
    TireConditions maxLongitudinalForces;

    // Calculate the limit for one tire.
    auto calculateFxMax = [&](double Fy, double Fz)
    {
        if (Fz <= 0.0 || longitudinalFriction <= 0.0 || lateralFriction <= 0.0)
            return 0.0;

        double maxFx = longitudinalFriction * Fz;
        double maxFy = lateralFriction * Fz;

        double lateralUtilization = std::abs(Fy) / maxFy;

        lateralUtilization = std::clamp(lateralUtilization, 0.0, 1.0);

        return maxFx * std::sqrt( 1.0 - lateralUtilization * lateralUtilization);
    };

    maxLongitudinalForces.frontLeft = calculateFxMax(lateralForces.frontLeft, verticalLoads.frontLeft);
    maxLongitudinalForces.frontRight = calculateFxMax(lateralForces.frontRight, verticalLoads.frontRight);
    maxLongitudinalForces.rearLeft = calculateFxMax(lateralForces.rearLeft, verticalLoads.rearLeft);
    maxLongitudinalForces.rearRight = calculateFxMax(lateralForces.rearRight, verticalLoads.rearRight);

    return maxLongitudinalForces;
}


VehicleAccelerations calculateVehicleAccelerations( const VehicleParameters& currentCar, const VehicleState& car, const TireConditions& lateralForces, double steeringAngle)
{
    VehicleAccelerations acceleration;

    // ==========================================
    // 1. Convert steering angle to radians
    // ==========================================

    double delta = steeringAngle * pi / 180.0;

    // ==========================================
    // 2. Calculate axle distances from CG
    // ==========================================

    double a = (currentCar.wheelbase / 2.0) - currentCar.gravityCenterOffset[0];
    double b = (currentCar.wheelbase / 2.0) + currentCar.gravityCenterOffset[0];

    // Half-track widths
    double frontHalfTrack = currentCar.frontTrackWidth / 2.0;
    double rearHalfTrack = currentCar.rearTrackWidth / 2.0;

    // ==========================================
    // 3. Transform tire forces into body axes
    // ==========================================

    // Front tires are steered by delta.
    // Longitudinal tire forces are currently zero.

    double frontLeftFx = -lateralForces.frontLeft * sin(delta);
    double frontRightFx = -lateralForces.frontRight * sin(delta);
    double frontLeftFy = lateralForces.frontLeft * cos(delta);
    double frontRightFy = lateralForces.frontRight * cos(delta);

    // Rear tires are not steered.
    double rearLeftFx = 0.0;
    double rearRightFx = 0.0;

    double rearLeftFy = lateralForces.rearLeft;
    double rearRightFy = lateralForces.rearRight;

    // ==========================================
    // 4. Calculate total vehicle forces
    // ==========================================

    double totalFx = frontLeftFx + frontRightFx + rearLeftFx + rearRightFx - car.dragForce;
    double totalFy = frontLeftFy + frontRightFy + rearLeftFy + rearRightFy;

    // ==========================================
    // 5. Calculate yaw moments about the CG
    // ==========================================

    // General equation:
    // Mz = x*Fy - y*Fx

    double frontLeftMoment = a * frontLeftFy - frontHalfTrack * frontLeftFx;
    double frontRightMoment = a * frontRightFy + frontHalfTrack * frontRightFx;

    double rearLeftMoment = -b * rearLeftFy - rearHalfTrack * rearLeftFx;
    double rearRightMoment = -b * rearRightFy + rearHalfTrack * rearRightFx;

    double yawMoment = frontLeftMoment + frontRightMoment + rearLeftMoment + rearRightMoment;

    // ==========================================
    // 6. Calculate vehicle accelerations
    // ==========================================

    // Body-fixed longitudinal acceleration
    acceleration.longitudinal = (totalFx / currentCar.mass) + car.yawRate * car.lateralVelocity;

    // Body-fixed lateral acceleration
    acceleration.lateral = (totalFy / currentCar.mass) - car.yawRate * car.longitudinalVelocity;

    // Angular acceleration about vertical axis
    acceleration.yaw = yawMoment / currentCar.yawInertia;

    return acceleration;
}


void updateVehicleMotion(VehicleState& car, const VehicleAccelerations& acceleration)
{
    car.longitudinalVelocity += acceleration.longitudinal * deltaTime;
    car.lateralVelocity += acceleration.lateral * deltaTime;
    double previousYawRate = car.yawRate;
    car.yawRate += acceleration.yaw * deltaTime;
    car.heading += 0.5 * (previousYawRate + car.yawRate) * deltaTime;
}

void aerodynamicEffect (const VehicleParameters& currentCar, VehicleState& car)
{
    car.dragForce = 0.5 * airDensity * currentCar.dragCoefficient * currentCar.frontalArea * pow(sqrt((pow(car.lateralVelocity , 2)) + (pow(car.longitudinalVelocity , 2))) , 2);
    car.downForce = 0.5 * airDensity * currentCar.downForceCoefficient * currentCar.liftSurfaceArea * pow(sqrt((pow(car.lateralVelocity , 2)) + (pow(car.longitudinalVelocity , 2))) , 2);
}

void calculateThrottle (DriverInput& control)    //Gradual keyboard throttle application
{
    const double throttleApplicationRate = 2.0;
    const double throttleReleaseRate = 3.0;
    //double throttle = 0;

    if (control.gasPressed)
        control.throttle += throttleApplicationRate * deltaTime;
    
    else
       control.throttle -= throttleReleaseRate * deltaTime; 
    
      control.throttle = std::clamp(control.throttle, 0.0, 1.0);
      
}
/*
void accelerate(DriverInput& control, VehicleState& car)
{
      const double accelerationRate = 5.0;
      car.longitudinalVelocity += accelerationRate * control.throttle * deltaTime;
}
*/

void applyDriverLongitudinalInput(const DriverInput& control, VehicleState& car)
{
    constexpr double maxAcceleration = 5.0; // m/s²
    constexpr double maxBraking = 8.0;      // m/s²
    double acceleration = maxAcceleration * control.throttle - maxBraking * control.brake;
    car.longitudinalVelocity += acceleration * deltaTime;
    // Temporary forward-only behavior.
    car.longitudinalVelocity = std::max(0.0, car.longitudinalVelocity);
}


void updatePosition(VehicleState& car)
{
    double worldVelocityX = car.longitudinalVelocity * cos(car.heading) - car.lateralVelocity * sin(car.heading);
    double worldVelocityY = car.longitudinalVelocity * sin(car.heading) + car.lateralVelocity * cos(car.heading);

    car.positionX += worldVelocityX * deltaTime;
    car.positionY += worldVelocityY * deltaTime;
}

void calculateBrake(DriverInput& control)
{
    constexpr double applicationRate = 4.0;
    constexpr double releaseRate = 6.0;

    if (control.brakePressed)
        control.brake += applicationRate * deltaTime;
    else
        control.brake -= releaseRate * deltaTime;

    control.brake =
        std::clamp(control.brake, 0.0, 1.0);
}

double steer(DriverInput& control,const VehicleParameters& currentCar)
{
    double steeringApplicationRate = 2.0;
    double steeringReleaseRate = 3.0;
    
    if (control.left && !control.right)
        control.steer += steeringApplicationRate * deltaTime; 
    
    else if (control.right && !control.left)
        control.steer -= steeringApplicationRate * deltaTime;
        
    else 
    {
        double releaseAmount = steeringReleaseRate * deltaTime;
        
        if (abs(control.steer) <= releaseAmount)
            control.steer = 0.0;
            
        else
            control.steer -= (control.steer / abs(control.steer))* releaseAmount;
    }
    
        control.steer = std::clamp(control.steer, -1.0, 1.0);
        double steeringAngle = control.steer * currentCar.steeringLockAngle;
        return steeringAngle;
}

//DriverInput control;

SimulationDiagnostics simulationStep(const VehicleParameters& currentCar, VehicleState& car, double steeringAngle, const CubicSpline& tireSpline, vector<vector<SplineData>>& tireData)
{
    SimulationDiagnostics diagnostics;

    // 1. Calculate aerodynamics
    aerodynamicEffect(currentCar, car);

    // 2. Calculate vertical loads
    diagnostics.verticalLoads = calculateVerticalForces(currentCar, car);

    // 3. Calculate slip angles
    diagnostics.slipAngles = SlipAngles(currentCar, car, steeringAngle);

    // 4. Calculate lateral forces
    diagnostics.lateralForces = calculateLateralForces(diagnostics.slipAngles, diagnostics.verticalLoads, tireSpline, tireData);

    // TEMPORARY DEBUGGING
    cout << "Slip FL: " << diagnostics.slipAngles.frontLeft << " | Force FL: " << diagnostics.lateralForces.frontLeft << '\n';

    // 5. Calculate vehicle accelerations
    diagnostics.accelerations = calculateVehicleAccelerations( currentCar, car, diagnostics.lateralForces, steeringAngle);
    
    diagnostics.lateralAcceleration = diagnostics.accelerations.lateral + car.yawRate * car.longitudinalVelocity;
    
    // 6. Update vehicle motion
    updateVehicleMotion(car, diagnostics.accelerations);

    // 7. Update vehicle position
    updatePosition(car);

    return diagnostics;
}

void processInput(DriverInput& control)
{
    control.gasPressed = sf::Keyboard::isKeyPressed(sf::Keyboard::Key::W);
    control.brakePressed = sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S);
    control.left = sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A);
    control.right = sf::Keyboard::isKeyPressed(sf::Keyboard::Key::D);
}


void drawVehicle(sf::RenderWindow& window, const VehicleState& car)
{
    sf::RectangleShape vehicle({60.f, 28.f});
    vehicle.setOrigin({30.f, 14.f});
    vehicle.setPosition({500.f, 350.f});
    vehicle.setFillColor(sf::Color(50, 150, 240));
    // World +Y points upward; screen +Y points downward.
    vehicle.setRotation(sf::degrees( static_cast<float>(-car.heading * 180.0 / pi)));
    window.draw(vehicle);
    // Small marker showing the front of the car.
    sf::CircleShape frontMarker(4.f);
    frontMarker.setOrigin({4.f, 4.f});
    frontMarker.setFillColor(sf::Color::Yellow);
    double angle = car.heading;
    frontMarker.setPosition({500.f + static_cast<float>(24.0 * cos(angle)), 350.f - static_cast<float>(24.0 * sin(angle))});
    window.draw(frontMarker);
}



int main()
{
    try
    {
        // ==========================================
        // 1. Initialize vehicle
        // ==========================================

        VehicleParameters currentCar;
        VehicleState car;
        DriverInput control;

        currentCar.mass = 1500.0;
        currentCar.wheelbase = 2.6;
        currentCar.frontTrackWidth = 1.6;
        currentCar.rearTrackWidth = 1.6;

        currentCar.gravityCenterOffset[0] = 0.0;
        currentCar.gravityCenterOffset[1] = 0.5;

        currentCar.steeringLockAngle = 30.0;
        currentCar.yawInertia = 2500.0;

        currentCar.dragCoefficient = 0.30;
        currentCar.downForceCoefficient = 0.10;
        currentCar.frontalArea = 2.2;
        currentCar.liftSurfaceArea = 2.2;
        currentCar.downForceFrontPercentage = 0.45;

        // Start moving so steering can be tested.
        car.longitudinalVelocity = 10.0;

        // ==========================================
        // 2. Load tire data
        // ==========================================

        auto tireData = dataFile("economy_tires");
        CubicSpline tireSpline(tireData);

        // ==========================================
        // 3. Create graphical window
        // ==========================================

        sf::RenderWindow window(sf::VideoMode({1000u, 700u}), "Vehicle Dynamics Simulator");
        window.setFramerateLimit(120);
        window.setKeyRepeatEnabled(false);

        sf::Clock clock;
        double accumulator = 0.0;
        bool paused = false;
        SimulationDiagnostics diagnostics;

        // ==========================================
        // 4. Real-time simulation loop
        // ==========================================

        while (window.isOpen())
        {
            // Measure elapsed real time.
            double frameTime = clock.restart().asSeconds();

            // Prevent enormous catch-up steps.
            frameTime = std::min(frameTime, 0.1);

            // ======================================
            // 5. Process window events
            // ======================================

            while (const auto event = window.pollEvent())
            {
                if (event->is<sf::Event::Closed>())
                    window.close();

                if (const auto* key = event->getIf<sf::Event::KeyPressed>())
                {
                    if (key->code == sf::Keyboard::Key::Escape)
                        window.close();

                    if (key->code == sf::Keyboard::Key::Space)
                        paused = !paused;

                    if (key->code == sf::Keyboard::Key::R)
                    {
                        car = VehicleState{};
                        control = DriverInput{};

                        car.longitudinalVelocity = 10.0;
                        diagnostics = SimulationDiagnostics{};
                        accumulator = 0.0;
                    }
                }
            }

            if (!window.isOpen())
                break;

            // ======================================
            // 6. Read continuous keyboard input
            // ======================================

            processInput(control);

            // ======================================
            // 7. Fixed-timestep physics
            // ======================================

            if (!paused)
            {
                accumulator += frameTime;

                while (accumulator >= deltaTime)
                {
                    calculateThrottle(control);
                    calculateBrake(control);
                    double steeringAngle = steer(control, currentCar);
                    // Temporary safety limit for
                    // the current tire model.
                    steeringAngle = std::clamp(steeringAngle, -3.0, 3.0);
                    // Temporary longitudinal model.
                    applyDriverLongitudinalInput(control, car);
                    // Limit speed during early testing.
                    car.longitudinalVelocity = std::min(car.longitudinalVelocity, 25.0);
                    diagnostics = simulationStep(currentCar, car, steeringAngle, tireSpline, tireData);
                    accumulator -= deltaTime;
                }
            }
            else
                accumulator = 0.0;

            // ======================================
            // 8. Calculate telemetry
            // ======================================

            double speed = std::hypot(car.longitudinalVelocity, car.lateralVelocity);

            double headingDegrees = car.heading * 180.0 / pi;

            double lateralG = diagnostics.lateralAcceleration / 9.81;

            // For now, display numerical telemetry
            // in the window's title bar.
            window.setTitle("Vehicle Dynamics | Speed: " + std::to_string(speed * 3.6) + " km/h | Heading: " + std::to_string(headingDegrees) + " deg | Yaw: " + std::to_string(car.yawRate) + " rad/s | Lateral G: " + std::to_string(lateralG) + (paused ? " | PAUSED" : ""));

            // ======================================
            // 9. Render the vehicle
            // ======================================

            window.clear(sf::Color(25, 28, 35));

            drawVehicle(window, car);

            window.display();
        }

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Simulation error: " << error.what() << '\n';
        return 1;
    }
}


