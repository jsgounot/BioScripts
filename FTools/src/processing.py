# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-02-17 13:45:38
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-02-17 14:01:07

from multiprocessing import Pool, cpu_count

class FunArgs() :
 
    def __init__(self, fun, * args, ** kwargs) :
        self.fun = fun
        self.args = args
        self.kwargs = kwargs
 
    def launch_fun(self) :
        return self.fun(* self.args, ** self.kwargs)
 
class Multiprocess() :
 
    @staticmethod
    def lambda_fun(farg) :
        return farg.launch_fun()
 
    def run(self, fargs, ncore=1) :
        if ncore == 1 : return (farg.launch_fun() for farg in fargs)
 
        pool = Pool(ncore)
        func = Multiprocess.lambda_fun
        fargs = ((farg, ) for farg in fargs)

        try :
            data = pool.starmap(func, fargs)
            pool.close()
            pool.join()
            return data
 
        except KeyboardInterrupt :
            print("Interrupt childs process")
            pool.terminate()
            sys.exit()
 
def mproc(fargs, ncore=1) :
    mp = Multiprocess()
    return mp.run(fargs, ncore=ncore)