The original source code is https://github.com/eic/EICrecon/tree/dirc_angleres/src/tests/pid_angleres by Matt posik

```bash
    #run reconstruction
    eicrecon \
      -Pplugins=pid_angleres \
      -PPidAngleRes:csv=${csv_out} \
      -PPidAngleRes:p=${p} \
      -PPidAngleRes:thetaL=${th_min} \
      -PPidAngleRes:thetaH=${th_max} \
      -Ppodio:output_file=${podio_out} \
      -Phistsfile=${hist_out} \
      "${in_root}"
```
