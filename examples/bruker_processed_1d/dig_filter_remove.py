import nmrglue as ng
import matplotlib.pyplot as plt
import os

data_dir  = os.path.join('bruker_data', '1')
procdata_dir = os.path.join(data_dir, 'pdata', '1')

dic, data = ng.bruker.read(data_dir)

#---time correction before FT
precorr_time =  ng.bruker.remove_digital_filter(dic, data)
precorr_frq = ng.proc_base.fft(precorr_time)
precorr_frq = ng.proc_autophase.autops(precorr_frq, 'peak_minima')

#---time correction after FT
postcorr_frq = ng.proc_base.fft(data)
postcorr_frq = ng.bruker.remove_digital_filter(dic, postcorr_frq, post_proc=True)
postcorr_frq = ng.proc_autophase.autops(postcorr_frq, 'acme')

#---data processed using TopSpin (3.5pl7)
bruker_frq = ng.bruker.read_pdata(procdata_dir)[1][::-1]

#---data processed using NMRPipe (preprocessed)
pipe_precorr_frq = ng.pipe.read('pre.ft')[1][::-1]
pipe_precorr_frq = ng.proc_autophase.autops(pipe_precorr_frq, 'peak_minima')

#---data processed using NMRPipe (post-procesing)
pipe_postcorr_frq = ng.pipe.read('post.ft')[1][::-1]
pipe_postcorr_frq = ng.proc_autophase.autops(pipe_postcorr_frq, 'acme') 

# note:
# 'acme' gave a better phase correction for the post-processing case
# 'peak_minima' gave a better phase correction for the pre-processing case

#---plot
fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(15, 5))

ax[0].plot(precorr_frq.real, label='NP='+str(precorr_frq.shape[-1]))
ax[0].set(title='nmrglue Pre-processing')

ax[1].plot(pipe_precorr_frq.real, label='NP='+str(pipe_precorr_frq.shape[-1]))
ax[1].set(title='NMRPipe Pre-processing')

ax[2].plot(pipe_postcorr_frq.real, label='NP='+str(pipe_postcorr_frq.shape[-1]))
ax[2].set(title='NMRPipe Post-processing')

ax[3].plot(bruker_frq.real, label='NP='+str(bruker_frq.shape[-1]))
ax[3].set(title='Processed using TopSpin')

ax[4].plot(postcorr_frq.real, label='NP='+str(postcorr_frq.shape[-1]))
ax[4].set(title='nmrglue Post-processing')


for axis in ax:
    axis.set(xlabel='Points', yticks=[])
    axis.legend()
    axis.invert_xaxis()

plt.tight_layout()
plt.savefig('dig_filter_remove.png')

