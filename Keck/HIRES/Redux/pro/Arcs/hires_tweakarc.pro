;+ 
; NAME:
; hires_tweakarc   
;     Version 1.1
;
; PURPOSE:
;  Allows user to edit the wavelength solution by hand
;
; CALLING SEQUENCE:
;  hires_tweakarc, twkfil, ordrs, templfil, QAFIL=, OSTR_FIL=
;   
;
; INPUTS:
;   twkfil   -  IDL save file containing the fit structure
;   ordrs    - Orders to tweak
;   templfil - Template file for archiving
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file is updated
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Summer-2005 Created by JXP 
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_tweakarc, twkfil, ordrs, templfil, _EXTRA=extra, $
                    QAFIL=qafil, OSTR_FIL=ostr_fil 


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_tweakarc, twkfil, ordrs, [templfil], QAFIL= [v1.0]'
      return
  endif 
  
;;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' 

  ;; Open line list
;  x_arclist, linlist, lines


  ;; Open templfil
  if keyword_set( TEMPLFIL ) then begin
      templfil = getenv('HIRES_CALIBS')+'/ARCS/'+templfil
      
;      tguess = guess_ordr
;      tfit = all_arcfit
;      tspec = sv_aspec
  endif

  x_tweakarc, twkfil, ordrs, templfil, QAFIL=qafil, OSTR_FIL=ostr_fil, $
    _EXTRA=extra

; All done
  print, 'hires_tweakarc: All done!'

  return
end
