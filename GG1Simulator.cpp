/*
 * File:   serverSim.cpp
 * Author: david tran
 * Version 1.2
 * Created on April 11, 2017, 11:42 AM
 * Purpose of this program is to simulate G/G/1 server queueing system.
 * The user will have the ability to decide what type of distribution for
 * interarrvial times and service times of jobs. The program shall output
 * statistics for the entire simulations depending on number of jobs the user
 * simulates. (minimum of 10000 jobs)
 * The program outputs identical statistics for the same number of jobs that have been served
 * this means that the job has not finish being served but another event/job has arrived. This
 * allows us to see how the queue changes.

 **User input variables**
 jMax = number of jobs for simulation
 iDA = iterarrival time distribution type
 aA,sA = range for interarrival time distribution
 iDS = for service time distribution type
 sA,sS = range for service time distribution
 **Queue variables**
 qA[] = Array holding arrival time of each job waiting in Queue
 qS[] = Array holding service time of jobs waiting in Queue
 in = next free place in the Queue
 out = job currently being served

 **Service and Timing variables**
 time = current time of simulation
 tNA = next arrival time
 tND = next departure time
 jobA = total jobs that have arrived
 jobS = total jobs that have been served
 tS = random service time
 tA = random interarrivaltime
 deltaT = time increment between events such as an arrival or departure
 nJob = number of jobs in entire service center

 **Statistics Variables**
 tBusy = total busy server time
 jobSec = job-seconds accumulated at the service center
 sDJob = standard deviation of number of jobs in system
 aTS = average service time during simulation
 sDTS = standard deviation of service time during simulation
 aTA = average interarrival time during simulation
 sDTA = standard deviation of interarrival time during simulation
 aRT = average response time during simulation
 sDRT = standard deviation of response time during simulation
 u = server utilization
 aJob = average number of jobs in queue during simulation
 y = throughput of jobs through server
 */
#include <cstdlib>
#include <iostream>
#include <math.h>
using namespace std;

//Function Prototypes
double normDisRNG(double rN, double a, double s),
uniformDisRNG(double rN, double a, double s),
uRN(),
expoDisRNG(double rN, double a),
arbDisRNG(double rN, double a, double s),
rng(int id, double rN, double a, double s);
int inc(int n);

//Const variable for queue size
int const qSize = 500;

//Main function of program
int main()
{
    //queue variables
    int in=0, out=0, jobA=0, jobS=0, nJob=0;
    int flag=1;
    double qA[qSize]={0.0}, qS[qSize]={0.0};

    //Statistics variables
    double tBusy=0.0, jobSec=0.0, sDJob=0.0, tNA=0.0, tND=0.0, aTS=0.0, sDTS=0.0, aTA=0.0, sDTA=0.0,
            aRT=0.0, sDRT=0.0, deltaT=0.0, time=0.0, tA=0.0, tS=0.0;

    //user input variables
    double aA=0.0, sA=0.0, aS=0.0, sS=0.0;
    int jMax, iDA, iDS;
    cout << "Directions for use of program G/G/1 simulator";
    cout << "jMax = number of jobs for simulation\n"
    << "iDA = iterarrival time distribution type please enter a number (1:Constant 2:Exponential 3:Uniform 4:Arbitrary 5:Normal)\n"
    << "aA,sA = range for interarrival time distribution\n"
    << "iDS = service time distribution type please enter a number (1:Constant 2:Exponential 3:Uniform 4:Arbitrary 5:Normal)\n"
    << "sA,sS = range for service time distribution\n\n";
    cout << "Enter inputs: JMAX IDA AA SA IDS AS SS" << endl;
    cin >> jMax >> iDA >> aA >> sA >> iDS >> aS >> sS;

    //Loop of main body
    while (flag)
    {

        //Event where time of next departure job is less than time of next arrival job, (a job will be served/depart)
        if (tND < tNA)
        {
          // calculate statistics variables for when a job has been served
            nJob = jobA - jobS;
            deltaT = tND - time;
            tBusy += deltaT;
            jobSec += (deltaT * nJob);
            sDJob += (deltaT * nJob) * nJob;
            time = tND;
            jobS = jobS + 1;
            aTS += (qS[out]);
            sDTS += (qS[out] * qS[out]);
            aRT += (time - qA[out]);
            sDRT = sDRT + (time - qA[out])*(time - qA[out]);

            //increment next index into queue
            out = inc(out);
            // checks to see if the queue is empty, if so tND will be calculated at next job arrival
            if (jobS == jobA)
                tND = tNA;
            else
                tND = time + qS[out];
            // checks if jobs served is = to maximum number of jobs to be simulated sets flag to 0 ends loop.
            if(jobS >= jMax)
                flag = 0;
        }
        // Event where a job has arrived to the service center
        else
        {
            nJob = jobA - jobS;
            if(nJob < qSize)
            {
                // Generate a random number for interarrival time for a job, loop until a non negative number is generated
                do
                {
                    tA = rng(iDA, tA, aA, sA);
                }
                while(tA <= 0);

                // generates a random number for time of service for a job, loop until a non negative number is generated
                do
                {
                    tS = rng(iDS, tS, aS, sS);
                }
                while(tS <= 0);

                //Calculate statistics variables for when a job has arrived to the service center
                jobA++;
                aTA += tA;
                sDTA += (tA*tA);

                //Checks to make sure that there are jobs in the service center to calculate these statistics
                if(nJob > 0)
                {
                    deltaT = tNA-time;
                    tBusy += deltaT;
                    jobSec += (deltaT*nJob);
                    sDJob += deltaT*(nJob*nJob);
                }
                else
                    tND = tNA + tS;

                //Calculate statistics variables for when a job has arrived to the service center continued
                time = tNA;
                tNA = time + tA;
                qA[in] = time;
                qS[in] = tS;
                //increment next index for job arrival
                in = inc(in);
            }
            else
                cout << "overflow error";
        }

        //Print all job statistics every 10000 jobs
        if (jobS%10000 == 0 && jobS > 0)
        {
        //Calculation to get statistics at certain time in the simulation
        double xaJob = jobSec/time;
        double u = tBusy/time;
        double xsDJob = sqrt(sDJob/(time)-(xaJob*xaJob));
        double xaTS = aTS/jobS;
        double xsDTS = sqrt(sDTS/(jobS) - (xaTS*xaTS));
        double xaTA = aTA/jobA;
        double xsDTA = sqrt(sDTA/(jobA) - (xaTA*xaTA));
        double xaRT = aRT/jobS;
        double xsDRT = sqrt(sDRT/(jobS) - (xaRT*xaRT));
        double y = jobS/time;
        double queueLength = u/(1-u);

        //Print statistics for current job divisible by 10000
        cout << "Jobs = " << jobS << endl << "U = " << u << endl << "AJob = " << xaJob << endl << "sDJob = " << xsDJob
        << endl << "ATS = " << xaTS << endl << "SDTS = " << xsDTS << endl << "aTA = "
        << xaTA << endl << "SDTA = " << xsDTA << endl << "ART = " << xaRT
        << endl << "SDRT = " << xsDRT << endl << "Y = " << y  << endl << "queueLength: " << queueLength << "\n\n\n";
        }
    }

    return 0;

}

/* Chooses the RNG distribution from user
 * @param id is the number used to select which RNG distribution to use
 *        a and s is the parameters for the distributions
 * @return a random number based on what function is chosen
 */
double rng(int id, double rN, double a, double s) {
    switch (id) // switch case to choose between different RNG distribution functions
    {
        case 5: // call to normal distribution RNG function
            return normDisRNG(rN, a, s);
        case 4: // call to arbitrary distribution RNG function
            return arbDisRNG(rN, a, s);
        case 3: // call to uniform distribution RNG function
            return uniformDisRNG(rN, a, s);
        case 2: // call to exponential distribution RNG function
            return expoDisRNG(rN, a);
        case 1: // in the case of constant value number
            return a;
        default:
            cout << "Please enter a correct distribution for RNG";
            exit(EXIT_FAILURE);
    }
}

//Function that returns a unique random number
double uRN() {
    return double(rand()) / RAND_MAX;
}

//Function returns a normally distributed random number
double normDisRNG(double rN, double a, double s) {
    double u = uRN();
    for (int i = 0; i < 10; i++) {
        u += uRN();
        rN = a + s * (u - 6.0);
    }
    return rN;
}

//Function returns uniformly distributed random number
double uniformDisRNG(double rN, double a, double s) {
    double u = uRN();
    rN = a + u * (s - a);
    return rN;
}

//Function returns exponentialy distributed random number
double expoDisRNG(double rN, double a) {
    double u = uRN();
    rN = -a * log(u);
    return rN;
}

//Function returns arbitrarily distributed random number
double arbDisRNG(double rN, double a, double s) {
    double u = uRN();
    rN = (sqrt(s * s + 2. * a * u) - s) / a;
    return rN;
}

// Function used to increment next index of queue for Ins and Outs of jobs
int inc(int n) {
    return ++n%qSize;
}
