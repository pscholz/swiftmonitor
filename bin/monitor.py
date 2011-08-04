#!/usr/bin/env python

from swiftmonitor import observation, config, mailer
import ftplib
import time
import numpy as np
import os
import re

def check_for_new():
  ftp = ftplib.FTP('heasarc.gsfc.nasa.gov')
  ftp.login()

  time_now = time.localtime()

  month_str = time.strftime('%Y_%m')  

  last_month_str = time.strftime('%Y_%m', time.strptime( '%s %s' % ( time_now.tm_year, time_now.tm_mon - 1 ), '%Y %m' ) )

  print 'Checking for new observations in /FTP/swift/data/obs/' + month_str +\
         ' and /FTP/swift/data/obs/' + last_month_str + ' ...'
  
  ftp.cwd('/FTP/swift/data/obs/' + month_str)
  obslist = ftp.nlst()
  ftp.cwd('/FTP/swift/data/obs/' + last_month_str)
  obslist.append( ftp.nlst() )

  targetids = np.loadtxt(config.monitor_code_path + 'targetids.list',dtype='S20')
  processed = np.loadtxt(config.monitor_code_path + 'processed.list',dtype='S20')
  obsids = []

  for obsid in obslist:
    if any( target in obsid for target in targetids ) and ( obsid not in processed ): 
      obsids.append(obsid)

  ftp.quit()

  return obsids

"""
def _download(ftp, obsids, month_str, last_month_str):

  wt_point_re = re.compile('*.xwtw2po_cl.evt*.') 

  for obsid in obsids:
    if obsid in ftp.nlst('/FTP/swift/data/obs/' + month_str):
      evt_lst = ftp.nlst('/FTP/swift/data/obs/' + month_str + '/' + obsid +\
                '/xrt/event/')
      for filn in evt_lst:
        if wt_point_re.match(filn):
"""       


def process(obslist):

  monitor_path = config.monitor_data_path
  if monitor_path[-1] != '/':
    monitor_path = monitor_path + '/'

  f = open(config.monitor_code_path + 'processed.list','w+')
  for obsid in obslist:

    try:

      obs_path = monitor_path + obsid
      os.mkdir(obs_path)
      ob = observation.Observation(obsid,obs_path)
      ob.download()
      ob.barycentre()
      ob.define_regions()
      ob.extract_region()
      ob.extract_spectrum()
      ob.correct_backscal()
      ob.save()

    except TableDownError as e:
      print e
      time.sleep(600)
      continue
    
    f.write(obsid + '\n')

  f.close() 


def main():
  obslist = check_for_new()

  if len(obslist) > 0:

    mail_str = 'New observations in archive:\n'
    for obsid in obslist:
      mail_str += '\t%s\n' % obsid
    mail_str += '\nWill now download and process them.'
    msg = mailer.MonitorMailer(mail_str, subject="New Observations")
    msg.send()  

    print mail_str

    process(obslist)

  else:
    print 'No new observations in archive.'

if __name__ == '__main__':
  main()
