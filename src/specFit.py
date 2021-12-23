#File for calculation and plot of fit

import numpy as np
from lmfit import Parameters, minimize, report_fit
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

#own modules
import dipole
import errorHandler as err
import util

#For tests
from pprint import pprint

class Fit:
  def __init__(self,config,ft,guess,Id,calcFlag='no'):
    #set id from input
    self.fitId = Id
    #read propagation time out of fourierTransform
    self.propTime = ft.propTime
    #read absolute error criterium out of config file
    self.fit_relerr_crit = config.fit_relerr_crit
    #read maximum of fit iterations out of config file
    self.fit_max_iter = config.fit_max_iter
    #read osci files out of guess
    self.osci = guess.osci
    #read guess out of guess-object
    self.guess = guess.guess
    #read fit_range
    self.fit_range = config.fit_range
    #calculate dw as the frequency step
    self.dw = (self.osci[0][len(self.osci[0][:,0])-1,0]-self.osci[0][0,0]) /\
              float((len(self.osci[0][:,0])-1))

    #run actuel fit
    self.makeFit(config,calcFlag) #saves fit results in self.fit_result

#-----------------------------------------------------------------------------#
#   Methods of the class Fit
#-----------------------------------------------------------------------------#
  def makeFit(self,config,calcFlag='no'):
    #Make an list for the guess, 
    # - if three files are fitted (calcFlag=='no'), then the list contains the
    #   guess for the x-, y- and z-direction
    # - if only one file is fitted (calcFlag=='trace'), then the list only contains
    #   one guess.
    
    # 1.) Make multifit or make single fit for trace
    #     The following routine call calculates the fit of the trace, as well!
    for i in range(self.fit_max_iter):
      if i+1 == self.fit_max_iter:
        print('Maximum number of iterations for fit exceeded for file ' + str(self.fitId+1))
      
      #Calculate the multifit (or in case of trace singlefit)
      self.fit_result = self.makeMultifit()

      #Calculate relative error between fit and raw data (in case of multifit the largest error is returned)
      self.fit_relerr = self.calcFitRelErr()
      print('fit_relerr = ' + str(self.fit_relerr) + ', fitId = ' + str(self.fitId))
      
      #stop if error criterium is fullfilled, else add a new line at maximum deviation
      #between fit and raw data
      if self.fit_relerr < self.fit_relerr_crit:
        print('Fit converged to file ' + str(self.fitId+1))
        break
      else:
        #add a new line in the spectrum and give a new guess at the position of the
        #maxiumum deviation between fit and raw data
        self.addNewLine()


    # 2.) Plot the result of the fit, the raw data and deviation
    self.plotFit(calcFlag)
    #TODO:
    # - Print fit results for every line in yaml file.

  #Routine for multifit
  def makeMultifit(self):
    #Make a list of w_i and f_i for all values of x-, y- and z-
    fit_params = Parameters()
    self.osci_range = np.abs(self.osci[0][0,0] - self.osci[0][len(self.osci[0][:,0])-1,0])
    for j in range(len(self.guess[0,1:])):
      off = int(j*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        fit_params.add('w_%i' % (i+1+off), value=self.guess[i,0])
        fit_params.add('f_%i' % (i+1+off), value=self.guess[i,j+1])
      
    #Make constraint that 'w_1' must equal 'w_(1+off)'
    for j in range(1,len(self.guess[0,1:])):
      off = int(j*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        fit_params['w_%i' % (i+1+off)].expr = 'w_%i' % (i+1)

    #calculate Fit
    #x is a list of numpy arrays with the energy axis of every direction
    #  in the case of 'trace' there is only x[0] as a numpy array of the energy
    x = []
    for j in range(len(self.guess[0,1:])):
      x.append(self.osci[j][:,0])
    #data is a list of the imaginary part of the fourier transform
    data = self.osci[0][:,2]
    for j in range(1,len(self.guess[0,1:])):
      data = np.append(data,self.osci[j][:,2])

    #Calculate fit
    out = minimize(self.objective, fit_params, args=(x,data))
    report_fit(out.params)

    #create fit_result as column_stack([w,f])
    fit_result = np.empty((len(self.guess[:,0]),0))
    w = np.zeros(len(self.guess[:,0]))
    for i in range(len(self.guess[:,0])):
      w[i] = out.params['w_%i' % (i+1)]
    fit_result = np.column_stack([fit_result,w])

    #number_direct is 1 for trace and 3 for multifit
    f = np.zeros(len(self.guess[:,0]))
    for j in range(len(self.guess[0,1:])):
      off = int(j*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        f[i] = out.params['f_%i' % (i+1+off)]

      fit_result = np.column_stack([fit_result,f])
        
    return fit_result

  #Routine for whole fit-function
  #x is a list of numpy arrays with the energy axis of every direction for the trace-fit
  #  there is only x[0]
  def fit_sinc(self, w, f, x):
    sinc = np.array([])
    for k in range(len(x)):
      s = np.zeros(len(x[k]))
      off = int(k*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        s += f[i+off]/w[i+off]*\
             np.sinc(self.propTime*(x[k]-w[i+off])/np.pi)*self.propTime
      
      sinc = np.append(sinc,s)
    
    return sinc


  #Routine for calculatin the total residual for fits of sinc-Function to osci_imag
  #x is an array appending self.osci[0][:,0], self.osci[1][:,0], self.osci[2][:,0]
  #data is an array appending self.osci[0][:,2], self.osci[1][:,2] and self.osci[2][:,2]. I. e. self.osci[:2][:,2]
  def objective(self, params, x, data):
    #read out of params w and f in arrays
    w = np.array([])
    f = np.array([])
    for j in range(len(self.guess[0,1:])):
      off = int(j*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        w = np.append(w,params['w_%i' % (i+1+off)])
        f = np.append(f,params['f_%i' % (i+1+off)])
    
    #make residual per data set
    resid = data - self.fit_sinc(w, f, x)

    #now flatten this to a 1D array, as minimize() needs
    return resid

  #Routine for calculation the relative error between fit and raw data (in case
  #of multifit the minimum error is returned)
  #absolute error is defined as err(dat,fit) = sqrt( [int(data-fit)^2 dw] / [fit-range] ) 
  #(integrals over fit range)
  def calcFitRelErr(self):
    rel_error = np.array([])
    for i in range(len(self.osci)):
      abs_err = np.sqrt(np.sum((self.osci[i][:,2] -\
                self.sinc(self.fit_result[:,0],self.fit_result[:,i+1],self.osci[i][:,0]))**2.)*\
                self.dw/np.abs(self.fit_range[1] - self.fit_range[0]))
      reference_err = np.sqrt(np.sum((self.osci[i][:,2])**2.)*\
                self.dw/np.abs(self.fit_range[1] - self.fit_range[0]))
      rel_error = np.append(rel_error,abs_err/reference_err)

    return np.amin(rel_error)

  #Routine for adding new line in guess for new fit at the maximum deviation between
  #fit and raw data
  def addNewLine(self):
    max_dev = np.array([])
    max_index = np.array([])
    for i in range(len(self.osci)):
      dev = self.osci[i][:,2] -\
            self.sinc(self.fit_result[:,0],self.fit_result[:,i+1],self.osci[i][:,0])
      max_dev = np.append(max_dev,np.amax(dev))
      max_index = np.append(max_index,np.argmax(dev))

    #Frequency of maximum deviation between fit and raw data
    w_max = self.osci[int(np.argmax(max_dev))][int(max_index[np.argmax(max_dev)]),0]
    
    #find corresponding amplitude(s) to the frequency
    f_max = np.array([])
    for i in range(len(self.osci)):
      abs_diff = np.abs(self.osci[i][:,0] - w_max)
      index_smallest_diff = abs_diff.argmin()
      f_max = np.append(f_max,self.osci[i][index_smallest_diff,2]*w_max/self.propTime)

    newLine = np.concatenate([[w_max],f_max])
    self.guess = np.vstack((self.guess,newLine))
    self.guess = self.guess[self.guess[:,0].argsort()]

#-----------------------------------------------------------------------------#
#   Methods for Making Plot
#-----------------------------------------------------------------------------#
  def plotFit(self,calcFlag):
    plt.ion()
    if calcFlag == 'no':
      fig = plt.figure()
      gs = fig.add_gridspec(len(self.osci), hspace=0)
      axs = gs.subplots(sharex=True, sharey=True)
      for i in range(len(self.osci)):
        axs[i].plot(self.osci[i][:,0],self.osci[i][:,2], label='Data in ' + util.getDir(i)
                    + '-direction')
        axs[i].legend(loc='best')
        axs[i].plot(self.osci[i][:,0],
                    self.sinc(self.fit_result[:,0],self.fit_result[:,i+1],self.osci[i][:,0]),
                    label='Fit in ' + util.getDir(i) + '-direction')
        axs[i].legend(loc='best')
        axs[i].plot(self.osci[i][:,0],
                    self.osci[i][:,2] -
                    self.sinc(self.fit_result[:,0],self.fit_result[:,i+1],self.osci[i][:,0]),
                    label='Error in ' + util.getDir(i) + '-direction')
        for j in range(len(self.fit_result[:,0])):
          #axs[i].arrow(self.fit_result[j,0], 0., 0., self.fit_result[j,i+1], length_includes_head=True,
          #             head_width=0.08, head_length=0.000002)
          axs[i].annotate("", xy=(self.fit_result[j,0],self.fit_result[j,i+1]*self.propTime/self.fit_result[j,0]),
                          xytext=(self.fit_result[j,0], 0.), arrowprops={'arrowstyle':'->'})
        axs[i].set_xlabel('Energy (Ry)')
        axs[i].legend(loc='best')
        axs[i].set_ylabel('$S(\hbar \omega)$')

      for ax in axs:
        ax.label_outer()

      fig.suptitle("Fit of dipole file " + str(self.fitId + 1))
      
      plt.show(block=False)

    elif calcFlag == 'trace':
      fig = plt.figure()
      plt.plot(self.osci[0][:,0],self.osci[0][:,2], label='Data trace')
      plt.plot(self.osci[0][:,0],
               self.sinc(self.fit_result[:,0],self.fit_result[:,1],self.osci[0][:,0]),
               label='Fit for trace')
      plt.plot(self.osci[0][:,0],
               self.osci[0][:,2] -
               self.sinc(self.fit_result[:,0],self.fit_result[:,1],self.osci[0][:,0]),
               label='Error for trace')
      for j in range(len(self.fit_result[:,0])):
          plt.annotate("", xy=(self.fit_result[j,0],self.fit_result[j,1]*self.propTime/self.fit_result[j,0]),
                       xytext=(self.fit_result[j,0], 0.), arrowprops={'arrowstyle':'->'})
      plt.legend(loc='best')
      plt.xlabel('Energy (Ry)')
      plt.ylabel('$S(\hbar \omega)$')
      fig.suptitle("Fit of trace")
      plt.show(block=False)
                    


  #Function for plotting the fit-results
  def sinc(self, w, f, x):
    s = np.zeros(len(x))
    for i in range(len(w)):
      s += f[i]/w[i]*np.sinc(self.propTime*(x-w[i])/np.pi)*self.propTime

    return s
      

