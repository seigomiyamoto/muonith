{
  "DETECTOR_PARAMETERS": {
    "name": "det_00"
  , "tf_read_bin_list": false
  // angle_unit should be Tangent, tangent, TANGENT,
  // or Degree, degree, DEGREE,
  // or Radian, radian, RADIAN
  , "angle_unit": "tangent"
  // if angle_unit = radian or degree,
  // 0 is north, 90 (PI/2) is east,
  // 180 (PI) is south, -90 (-PI/2) is west
  // and should be 90 > tymax > tymin >= 0
  , "nbinx": $nbinx
  , "txmin": $txmin
  , "txmax": $txmax
  , "nbiny": $nbiny
  , "tymin": $tymin
  , "tymax": $tymax
  , "length_hori": $length_hori
  , "length_vert": $length_vert
  , "length_dept": $length_dept
  , "n_unit": $n_unit
  // rotation type LOCAL or GLOBAL
  , "rotation_type": "LOCAL"
  , "yaw_deg": -31.75
  , "roll_deg": 0.0
  , "pitch_deg": 0.0
  , "x": 386.4
  , "y": -327.9
  , "z": 335.6
  , "days": $days
  , "n_reserve_vec_tf_in_PL": 20
  , "path_eff_table" : "../detparams/$spec_label/eff_table_sample00.tmp"
  // Analytic efficiency model json5; takes precedence over path_eff_table
  // when set (e.g. "../detparams/eff_model/small_uncertainty.json5").
  , "path_eff_model" : "none"
  }
}
