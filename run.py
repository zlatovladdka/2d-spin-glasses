#!/usr/bin/python3

import datetime
import numpy as np
import subprocess

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
sim_group.add_argument("--seed", help="Random seed for J-lattice generation", type=int)
sim_group.add_argument("-d", "--dist", help="Distribution for J-lattice generation: 0 - Normal, 1 - uniform plus-minus J", type=int, default=0)
mc_group = parser.add_argument_group("Monte-Carlo")
mc_group.add_argument("--Ntherm", help="Number of steps for thermalization", type=int, default=10000)
mc_group.add_argument("--Nmc", help="Number of steps used for measurements", type=int, default=10000)
corr_group = parser.add_argument_group("Autocorrelator")
corr_group.add_argument("--tmax", help="Maximal time to measure autocorrelation function", type=float, default=5.0)
corr_group.add_argument("--dt", help="Time step for autocorrelation function", type=float, default=0.5)
corr_group.add_argument("--bsize", help="Bin size for averaging over spin flips", type=int, default=0)
db_group = parser.add_argument_group("Database")
db_group.add_argument("-H", "--host", help="MongoDB hostname", type=str, default='c1.itp.ac.ru')
db_group.add_argument("-P", "--port", help="MongoDB port", type=int, default=27017)
prl_group = parser.add_argument_group("Parallelization")
prl_group.add_argument("--Nproc", help="Number of parallel processes", type=int, default=1)

args = parser.parse_args()

Nx = args.size[0]
Ny = args.size[1]
temp = args.temperature
rand_seed = args.seed
tp = args.type
distType = args.dist

Ntherm = args.Ntherm
Nmc = args.Nmc

tmax = args.tmax
dt = args.dt
bin_size = args.bsize

n = args.Nproc

options = ["/home-parma/vtemkin/code/2d-spin-glasses/mc", "-x", str(Nx), "-y", str(Ny), "--temp", str(temp), "--dist", str(distType), "--therm", str(Ntherm), "--time", str(Nmc), "--autocorr", str(tmax), "--autocorr-dt", str(dt), "--bin-size", str(bin_size)]

flag = False
if str(tp).lower() == "random":
    options += ["-g"]
    flag = True


def run_mc(opts):
    if args.seed is not None:
    	rand_seed = args.seed
    else:
	    rand_seed = int(datetime.datetime.now().timestamp() * 1e6)
    if "--seed" not in opts:
        opts += ["--seed", str(rand_seed)]
    else:
        opts[-1] = str(rand_seed)

    client = pymongo.MongoClient(args.host, args.port)
    data = subprocess.check_output(opts, stderr=subprocess.DEVNULL)
    data = list(map(float, data.split()))
    print("Got system (E={0:.3f}, M={1:.3f})".format(data[0], data[2]))
    result = {"Nx": Nx,
              "Ny": Ny,
              "T": temp,
              "Seed": rand_seed,
              "Distribution": distType,
              "Ntherm": Ntherm,
              "Nmc": Nmc,
              "Bin_size": bin_size,
              "Energy": data[0],
              "Energy2": data[1],
              "Magnetization": data[2],
              "Magnetization2": data[3],
              "corr_times": data[4::2],
              "corr_values": data[5::2]}
    print(data[4::2])
    print(data[5::2])
    while True:
        if flag:
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

while True:
    run_mc(options)

