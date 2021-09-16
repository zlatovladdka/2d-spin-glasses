#!/usr/bin/python3

import numpy as np
import subprocess
from multiprocessing import Pool, Process

import pymongo
import gridfs
from bson.binary import Binary
from pickle import dumps, loads

compress = lambda x: Binary(dumps(x, protocol = 2), subtype = 128)

### Loading parameters
import argparse
parser = argparse.ArgumentParser(description="Ising model solver")
sim_group = parser.add_argument_group("Simulation")
sim_group.add_argument("-s", "--size", help="Size of the system", type=int, nargs=2, required=True)
sim_group.add_argument("-T", "--temperature", help="Temperature", type=float, required=True)
sim_group.add_argument("--type", help="Type of couplings: const or random", type=str, required=True)
sim_group.add_argument("--seed", help="Random seed for J-lattice generation", type=int, required=True)
mc_group = parser.add_argument_group("Monte-Carlo")
mc_group.add_argument("--Ntherm", help="Number of steps for thermalization", type=int, default=10000)
mc_group.add_argument("--Nmc", help="Number of steps used for measurements", type=int, default=10000)
corr_group = parser.add_argument_group("Autocorrelator")
corr_group.add_argument("--tmax", help="Maximal time to measure autocorrelation function", type=float, default=5.0)
corr_group.add_argument("--dt", help="Time step for autocorrelation function", default=0.5)
db_group = parser.add_argument_group("Database")
db_group.add_argument("-H", "--host", help="MongoDB hostname", type=str, default='c1.itp.ac.ru')
db_group.add_argument("-P", "--port", help="MongoDB port", type=int, default=27017)

args = parser.parse_args()

Nx = args.size[0]
Ny = args.size[1]
temp = args.temperature
rand_seed = args.seed
type = args.type
if str(type) == "random":
    flag = 1
else:
    flag = 0

Ntherm = args.Ntherm
Nmc = args.Nmc

tmax = args.tmax
dt = args.dt

opts = ["./new_ising.exe", "-x", str(Nx), "-y", str(Ny), "--temp", str(temp), "--therm", str(Ntherm), "--time", str(Nmc), 
        "--autocorr", str(tmax), "--autocorr-dt", str(dt), "-g" * flag + " ", "--seed", str(rand_seed)]

client = pymongo.MongoClient(args.host, args.port)

def run_mc(seed):
    data = subprocess.check_output(opts, stderr=subprocess.DEVNULL)
    data = list(map(float, data.split()))
    print("Got system (E={0:.3f}, M={1:.3f})".format(data[0], data[2]))
    result = {"Nx": Nx,
              "Ny": Ny,
              "T": temp,
              "Seed": seed,
              "Ntherm": Ntherm,
              "Nmc": Nmc,
              "Energy": data[0],
              "Energy2": data[1],
              "Magnetization": data[2],
              "Magnetization2": data[3],
              "corr_times": data[4::2],
              "corr_values": data[5::2]}
    print(data[4::2])
    print(data[5::2])
    while True:
        if flag == 1:
            try:
                client.numerics.glass.insert_one(result)
            except (pymongo.errors.ConnectionFailure, pymongo.errors.ServerSelectionTimeoutError) as e:
                print("Connection to DB failed, retrying...")
            else:
                print("Data inserted successfully to glass table!")
                break
        else:
            try:
                client.numerics.ising.insert_one(result)
            except (pymongo.errors.ConnectionFailure, pymongo.errors.ServerSelectionTimeoutError) as e:
                print("Connection to DB failed, retrying...")
            else:
                print("Data inserted successfully to ising table!")
                break

if __name__ == "__main__":
    while True:
        with Pool(6) as p:
            p.map(run_mc, [rand_seed]*6)