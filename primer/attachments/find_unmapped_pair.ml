#!/usr/bin/env ocamlscript
Ocaml.packs := [ "pcre"; "unix"; "extlib"; "str" ];;
Ocaml.ocamlflags := [ "-g" ];;
--
open Printf
open Arg

module Strand =
struct
  type t = PLUS | MINUS
  let of_string s = 
    if s = "+" then PLUS
    else if s = "-" then MINUS
    else raise (Invalid_argument "Strand.of_string")
end

module Hit = struct
  type t = {line:string; id:string; read_id:int; chr:string; pos:int; strand:Strand.t; ins:int; flag:int; seq:string; qual:string}

  let of_line line =
    match Pcre.split line with
    | s::chr::pos::strand::ins::flag::
      map_qual::se_qual::alt_qual::mismatches::sum_mis_qual::m0::m1::read_length::
      seq::qual::[] -> begin
        ignore [map_qual;se_qual;alt_qual;mismatches;sum_mis_qual;m0;m1;read_length;];
        match Pcre.split ~pat:"/" s with
        | id::read_id::[] -> begin
          {line=""; id=id; read_id=int_of_string read_id; chr=chr; 
           pos=int_of_string pos; strand=Strand.of_string strand;
           ins=int_of_string ins; flag=int_of_string flag;
           seq=seq; qual=qual}
          end
        | _ -> printf "failure parsing header: Hit.of_line %s\n" s; assert false
      end
    | _ -> printf "failure parsing line: Hit.of_line %s\n" line; assert false
    
  let has_same_key t1 t2 = t1.id = t2.id

  let pair_on_same_chromosome t1 t2 = t1.chr = t2.chr

  let pair_on_same_strand t1 t2 =
    (t1.strand = Strand.PLUS && t2.strand = Strand.PLUS) ||
    (t1.strand = Strand.MINUS && t2.strand = Strand.MINUS)

  let pair_reversed t1 t2 =
    (t1.strand = Strand.PLUS && t2.strand = Strand.MINUS &&
     t1.pos > t2.pos) ||
    (t2.strand = Strand.PLUS && t1.strand = Strand.MINUS &&
     t2.pos > t1.pos)

  let distance t1 t2 = abs (t1.pos - t2.pos)

  let is_unmapped t = (t.flag land 192) != 0

  let mate_is_unmapped t = (t.flag land 64) != 0

  let mate_mapped_to_different_chromosome t = (t.flag land 32) != 0

  let to_fastq_string t = 
    sprintf "@%s/%d\n%s\n+\n%s\n" t.id t.read_id t.seq t.qual
end

(* ****** *)

let last_hit = ref None
let lower_size_limit = ref 10
let upper_size_limit = ref 300
let output_circos = ref false
let output_table = ref false
let output_unmapped = ref false
(* let circos_file = ref "" *)
let files = ref []
let gff = ref ""

let flags = align
[
("-l", Set_int lower_size_limit, sprintf "INT lower size limit (default = %d)" !lower_size_limit);
("-L", Set_int upper_size_limit, sprintf "INT upper size limit (default = %d)" !upper_size_limit);
("-circos", Set output_circos, sprintf "output circos link file (default: %s)" (string_of_bool !output_circos));
("-table", Set output_table, sprintf "output a table (default: %s)" (string_of_bool !output_table));
("-unmapped", Set output_unmapped, sprintf "output unmapped reads and their mates (default: %s)" (string_of_bool !output_unmapped));
("-gff_filter", Set_string gff,  
             "file a file in GFF format specifying genomic features and their genomic 
                    coordinates; reads that are mapped to the specified regions will 
                    be printed.  For paired-end reads, only one of the mate pairs 
                    need to be mapped");
]

let anon_func x = files := x::(!files)

let () = parse flags anon_func ""

let default_print ~last_hit cur_hit line =
  if cur_hit.Hit.strand = Strand.MINUS then begin
    print_string ("-\t" ^ last_hit.Hit.line ^ "\n");
    print_string ("-\t" ^ line ^ "\n");
    flush_all  ()
  end else begin
    print_string ("+\t" ^ line ^ "\n");
    print_string ("+\t" ^ last_hit.Hit.line ^ "\n");
    flush_all  ()
  end

(* TODO -- should handle the situation where two mate pairs are mapped to different chromosomes *)
let print_circos ~last_hit cur_hit _ =
  printf "TODO -- should handle the situation where two mate pairs are mapped to different chromosomes\n";
  let (* chr1 = last_hit.Hit.chr  and
      chr2 =  cur_hit.Hit.chr  and *)
      s1 = last_hit.Hit.strand and
      s2 =  cur_hit.Hit.strand and
      p1 = last_hit.Hit.pos    and
      p2 =  cur_hit.Hit.pos    in
  let style = 
    if s1 = s2 then 
      if s1 = Strand.PLUS then "\tcolor=red" 
      else "\tcolor=blue" 
    else 
      if (s1 = Strand.PLUS && p1 > p2) ||
         (s2 = Strand.PLUS && p2 > p1) then
        "\tcolor=green"
      else
        if Hit.distance last_hit cur_hit < !lower_size_limit then
          "\tcolor=orange"
        else ""
  in
  printf "%s\t%s\t%d\t%d%s\n" cur_hit.Hit.id cur_hit.Hit.chr cur_hit.Hit.pos cur_hit.Hit.pos style;
  printf "%s\t%s\t%d\t%d%s\n" last_hit.Hit.id last_hit.Hit.chr last_hit.Hit.pos last_hit.Hit.pos style

let print_table ~last_hit cur_hit _ =
  let chr1 = last_hit.Hit.chr  and
      chr2 =  cur_hit.Hit.chr  and
      s1 = last_hit.Hit.strand and
      s2 =  cur_hit.Hit.strand and
      p1 = last_hit.Hit.pos    and
      p2 =  cur_hit.Hit.pos    in
  let style = 
    if chr1 <> chr2 then "5"
    else
    begin
      if s1 = s2 then 
       begin
        if s1 = Strand.PLUS then "3"
        else "4"
       end
      else (* opposite strand *)
       begin
        if (s1 = Strand.PLUS && p1 > p2) ||
           (s2 = Strand.PLUS && p2 > p1) then
          "2"
        else "1"
       end
    end
  in
  let p1,p2,seq1,seq2,chr1,chr2 =
    match s1, s2 with
    | Strand.PLUS , Strand.PLUS  -> 
        if p1 < p2 then 
             (p1,p2,last_hit.Hit.seq,cur_hit.Hit.seq,chr1,chr2) 
        else (p2,p1,cur_hit.Hit.seq,last_hit.Hit.seq,chr2,chr1)
    | Strand.MINUS, Strand.MINUS -> 
        if p1 > p2 then 
             (p1,p2,last_hit.Hit.seq,cur_hit.Hit.seq,chr1,chr2) 
        else (p2,p1,cur_hit.Hit.seq,last_hit.Hit.seq,chr2,chr1)
    | Strand.PLUS , Strand.MINUS -> 
        (p1,p2,last_hit.Hit.seq,cur_hit.Hit.seq,chr1,chr2) 
    | Strand.MINUS, Strand.PLUS  ->
        (p2,p1,cur_hit.Hit.seq,last_hit.Hit.seq,chr2,chr1)
  in
  printf "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n" cur_hit.Hit.id p1 p2 style seq1 seq2 chr1 chr2;
  flush_all ()

let print_unmapped ~last_hit cur_hit _ =
    print_string (Hit.to_fastq_string last_hit);
    print_string (Hit.to_fastq_string cur_hit);
    flush_all ()

let f_print =
  if !output_table then print_table
  else if !output_circos then print_circos
  else default_print

(* *** *)
let regions = 
  if !gff = "" then [||]
  else
    let regions = DynArray.create () in
    let fin = open_in !gff in
    let () = 
      try
        while true do
          let line = input_line fin in
          let fld = Array.of_list (Str.split (Str.regexp "[\t\n]") line) in
          if Array.length fld < 5 then ()
          else
            if String.get fld.(0) 0 = '#' then () 
            else
              let chr = fld.(0) and pos1 = int_of_string fld.(3) and pos2 = int_of_string fld.(4) in
              DynArray.add regions (chr,pos1,pos2)
        done
      with End_of_file -> ()
    in
    let regions = DynArray.to_array regions in
    let () =  Array.sort (fun (a,b,_) (a',b',_) -> 
                 match compare a a' with
                 | 0 -> compare b b' 
                 | x -> x ) regions
    in regions
  
let rec _find regions chr pos i j =
  if i >= j then raise Not_found
  else
    let k = (i+j)/2 in
    let a,b,c = regions.(k) in
    match compare chr a, pos >= b, pos <= c with
    | 0, true, true ->  a,b,c
    | -1, _, _ -> _find regions chr pos 0 (k-1)
    | 0, false, _ -> _find regions chr pos 0 (k-1)
    | 1,  _, _ -> _find regions chr pos (k+1) j
    | 0,  _, false -> _find regions chr pos (k+1) j
    | _, _, _ -> raise Not_found

let find regions chr pos = _find regions chr pos 0 (Array.length regions)

let skip_check = !gff = ""

let check_filter hit1 hit2 =
  let chr1 = hit1.Hit.chr and
      pos1 = hit1.Hit.pos and
      chr2 = hit2.Hit.chr and
      pos2 = hit2.Hit.pos 
  in
  try begin
    ignore (find regions chr1 pos1); true
  end
  with Not_found -> 
    try
      begin
        ignore (find regions chr2 pos2); true
      end
    with Not_found -> false


let line_number = ref 0

let do_line line =
  let cur_hit = Hit.of_line line in
  let () = begin
    match !last_hit with
    | Some last_hit -> begin
        if (Hit.has_same_key last_hit cur_hit) then 
        begin
          if !output_unmapped then (* find unmapped (singly mapped) pairs in the original maq map *)
          (
            if (skip_check || check_filter last_hit cur_hit) && (Hit.is_unmapped last_hit || Hit.is_unmapped cur_hit) then 
              print_unmapped ~last_hit cur_hit line
            else ()
          )
          else (* find anomalous paired reads *)
          (
            if Hit.is_unmapped last_hit || Hit.is_unmapped cur_hit then 
              ()
            else if (skip_check || check_filter last_hit cur_hit) && 
               ((not (Hit.pair_on_same_chromosome last_hit cur_hit)) || (Hit.pair_reversed last_hit cur_hit) || (Hit.pair_on_same_strand last_hit cur_hit) ||
                 (Hit.distance last_hit cur_hit > !upper_size_limit) ||
                 (Hit.distance last_hit cur_hit < !lower_size_limit)) then
              f_print ~last_hit cur_hit line
            else ()
          )
        end 
        else (* don't have the same key as the previous hit *) 
          ()
      end
    | _ -> ()
  end
  in last_hit := Some cur_hit
  
(*
let file = "s_2@1.unmapped.map25"
let cmd = "maq mapview " ^ file ^ " | sort "

let fin = Unix.open_process_in cmd
*)

let () =
  List.iter
    (fun file ->
      (* e.g. "s_2@1.unmapped.map25" *)
      let sorted_file = file ^ ".sorted" and 
          tmp_sorted_file = file ^ ".sorted.tmp" in
      let cmd = "if [ -e "^sorted_file^" ]; then cat "^sorted_file^"; else maq mapview "^file^" | sort >"^tmp_sorted_file^"; mv "^tmp_sorted_file^" "^sorted_file^"; cat "^sorted_file^"; fi" in 
      let fin = Unix.open_process_in cmd in
      let () = Pcre.foreach_line ~ic:fin do_line in
        close_in fin 
    ) !files

let test () =
  let l = "HWI-EAS00184:2:2:100:1461#0/1	fc40	207316	+	193	18	0	0	0	0	0	2	0	36	TTGAGGCCAACGCCCATAATGCGGGCGGtTgccCgG	BBBBA?BBBBABBBAABBBAA??=?AAB4>92;=6=" in
  let () = do_line l in
  let l = "HWI-EAS00184:2:2:100:1461#0/2	fc40	207473	-	-193	18	0	0	0	0	0	2	0	36	CGGCAGTGCTTTTGCCGTTACGCACCACCCCGTcaG	BABAB>A@@BBBBABBABBBAAB=BBAAABBBA;.<" in
  let () = do_line l in
    ()

