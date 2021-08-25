#!/usr/bin/python3

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

Ntherm = args.Ntherm
Nmc = args.Nmc

tmax = args.tmax
dt = args.dt

opts = ["./mc", str(Nx), str(Ny), str(temp), str(Ntherm), str(Nmc), str(tmax), str(dt)]

client = pymongo.MongoClient(args.host, args.port)

while True:
    data = subprocess.check_output(opts, stderr=subprocess.DEVNULL)
    print("Got system (E={0:.3f}, M={1:.3f})".format(data[0], data[2]))
    data = list(map(float, data.split()))
    result = {"Nx": Nx,
              "Ny": Ny,
              "T": temp,
              "Ntherm": Ntherm,
              "Nmc": Nmc,
              "Energy": data[0],
              "Energy2": data[1],
              "Magnetization": data[2],
              "Magnetization2": data[3],
              "corr_times": data[4::2],
              "corr_values": data[5::2]}
    while True:
        try:
            client.numerics.ising.insert_one(result)
        except (pymongo.errors.ConnectionFailure, pymongo.errors.ServerSelectionTimeoutError) as e:
            print("Connection to DB failed, retrying...")
        else:
            print("Data inserted successfully!")
            break
