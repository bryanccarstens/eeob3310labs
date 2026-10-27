text_color <- function(hex) {
  r <- strtoi(substr(hex, 2, 3), 16)
  g <- strtoi(substr(hex, 4, 5), 16)
  b <- strtoi(substr(hex, 6, 7), 16)
  lum <- (0.299 * r + 0.587 * g + 0.114 * b) / 255
  ifelse(lum > 0.5, "#000000", "#ffffff")
}

aa_colors <- c(
  A = "#414357", R = "#3c4e61", N = "#386269", D = "#2b8273",
  C = "#4da361", E = "#9bcf70", Q = "#c8e687", G = "#fcffa3",
  H = "#ffe375", I = "#ffbe5c", L = "#ff954f", K = "#ff643d",
  M = "#f52c2c", F = "#b83359", P = "#852970", S = "#59226b",
  T = "#3d2352", W = "#4a2d75", Y = "#503785", V = "#566db0",
  "-" = "#1a1a2e", X = "#333333"
)

nt_colors <- c(
  A = "#4da361", T = "#f52c2c", G = "#ffbe5c",
  C = "#566db0", N = "#333333", "-" = "#1a1a2e"
)

compute_stats <- function(seqs) {
  mat <- do.call(rbind, strsplit(seqs, ""))
  ncols <- ncol(mat)
  identical <- sum(apply(mat, 2, function(col) !any(col == "-") && length(unique(col)) == 1))
  mismatched <- sum(apply(mat, 2, function(col) !any(col == "-") && length(unique(col)) > 1))
  gapped <- sum(apply(mat, 2, function(col) any(col == "-")))
  list(
    identical = round(100 * identical / ncols, 1),
    mismatched = round(100 * mismatched / ncols, 1),
    gapped = round(100 * gapped / ncols, 1)
  )
}

build_alignment_html <- function(seqs, uid, colors, label_width = 210, img_uris = NULL) {
  names_list <- names(seqs)
  stats <- compute_stats(seqs)
  nseqs <- length(seqs)
  seq_len <- nchar(seqs[[1]])
  cell_w <- 18
  cell_h <- 28
  canvas_h <- nseqs * cell_h
  max_cols_per_tile <- 1500
  n_tiles <- ceiling(seq_len / max_cols_per_tile)
  
  color_json <- paste0("{", paste(sprintf('"%s":"%s"', names(colors), colors), collapse = ","), "}")
  fg_json <- paste0("{", paste(sprintf('"%s":"%s"', names(colors), sapply(colors, text_color)), collapse = ","), "}")
  seq_json <- paste0("[", paste(sprintf('"%s"', seqs), collapse = ","), "]")
  
  label_html <- ""
  for (i in seq_along(seqs)) {
    sp <- names_list[[i]]
    img_tag <- ""
    if (!is.null(img_uris) && sp %in% names(img_uris) && nchar(img_uris[[sp]]) > 0) {
      img_tag <- paste0('<img src="', img_uris[[sp]], '" style="width:24px; height:24px; object-fit:contain;">')
    }
    label_html <- paste0(label_html,
                         '<div style="display:flex; align-items:center; justify-content:flex-end; height:', cell_h, 'px; padding-right:8px; gap:10px; cursor:pointer;" onclick="toggleHL_', uid, '(', i - 1, ')" data-idx="', i - 1, '">',
                         img_tag,
                         '<span style="font-weight:bold; color:#333; font-size:13px; font-family:monospace; white-space:nowrap;">', sp, '</span>',
                         '</div>')
  }
  
  canvas_tags <- ""
  for (t in seq_len(n_tiles)) {
    col_start <- (t - 1) * max_cols_per_tile
    col_end <- min(t * max_cols_per_tile, seq_len)
    tile_w <- (col_end - col_start) * cell_w
    canvas_tags <- paste0(canvas_tags,
                          '<canvas id="cv-', uid, '-', t, '" width="', tile_w, '" height="', canvas_h,
                          '" style="display:inline-block; vertical-align:top;"></canvas>')
  }
  
  paste0(
    '<div style="display:flex; gap:24px; padding:12px 0; font-family:monospace; font-size:14px; margin-bottom:8px;">',
    '<div style="display:flex; align-items:center; gap:8px;"><span style="width:10px; height:10px; border-radius:50%; background:#4da361; display:inline-block;"></span>',
    '<span style="color:#666;">Identical:</span><span style="color:#333; font-weight:bold;">', stats$identical, '%</span></div>',
    '<div style="display:flex; align-items:center; gap:8px;"><span style="width:10px; height:10px; border-radius:50%; background:#ff643d; display:inline-block;"></span>',
    '<span style="color:#666;">Mismatched:</span><span style="color:#333; font-weight:bold;">', stats$mismatched, '%</span></div>',
    '<div style="display:flex; align-items:center; gap:8px;"><span style="width:10px; height:10px; border-radius:50%; background:#ccc; border:1px solid #999; display:inline-block;"></span>',
    '<span style="color:#666;">Gapped:</span><span style="color:#333; font-weight:bold;">', stats$gapped, '%</span></div>',
    '</div>',
    
    '<div style="display:flex; position:relative;" id="aln-wrap-', uid, '">',
    
    '<div style="flex-shrink:0; width:', label_width, 'px; background:#fff; z-index:2;" id="label-col-', uid, '">',
    label_html,
    '</div>',
    
    '<div style="flex:1; min-width:0; overflow-x:auto; overflow-y:hidden;" id="scroll-', uid, '">',
    '<div style="white-space:nowrap; font-size:0; line-height:0;">',
    canvas_tags,
    '</div>',
    '</div>',
    
    '</div>',
    
    '<script>',
    '(function(){',
    'var uid="', uid, '";',
    'var cw=', cell_w, ',ch=', cell_h, ',totalCols=', seq_len, ',nseqs=', nseqs, ';',
    'var maxCPT=', max_cols_per_tile, ',nTiles=', n_tiles, ';',
    'var seqs=', seq_json, ';',
    'var bg=', color_json, ';',
    'var fg=', fg_json, ';',
    'var hl=new Set();',
    
    'function drawTile(t){',
    '  var cv=document.getElementById("cv-"+uid+"-"+t);',
    '  if(!cv)return;',
    '  var ctx=cv.getContext("2d");',
    '  var colStart=(t-1)*maxCPT;',
    '  var colEnd=Math.min(t*maxCPT,totalCols);',
    '  ctx.clearRect(0,0,cv.width,cv.height);',
    '  ctx.textAlign="center";ctx.textBaseline="middle";',
    '  ctx.font="bold 12px monospace";',
    '  var anyHL=hl.size>0;',
    '  for(var i=0;i<nseqs;i++){',
    '    var dimmed=anyHL&&!hl.has(i);',
    '    ctx.globalAlpha=dimmed?0.25:1;',
    '    var s=seqs[i];',
    '    for(var j=colStart;j<colEnd;j++){',
    '      var x=(j-colStart)*cw;',
    '      var c=s.charAt(j).toUpperCase();',
    '      ctx.fillStyle=bg[c]||"#333";',
    '      ctx.fillRect(x,i*ch,cw,ch);',
    '      ctx.fillStyle=fg[c]||"#fff";',
    '      ctx.fillText(c,x+cw/2,i*ch+ch/2+1);',
    '    }',
    '  }',
    '  ctx.globalAlpha=1;',
    '}',
    
    'function drawAll(){for(var t=1;t<=nTiles;t++)drawTile(t);}',
    
    'window.toggleHL_', uid, '=function(idx){',
    '  if(hl.has(idx))hl.delete(idx);else hl.add(idx);',
    '  drawAll();',
    '  var labels=document.querySelectorAll("#label-col-"+uid+" [data-idx]");',
    '  labels.forEach(function(l){',
    '    var li=parseInt(l.getAttribute("data-idx"));',
    '    l.style.opacity=(hl.size>0&&!hl.has(li))?"0.3":"1";',
    '  });',
    '};',
    
    'if(document.readyState==="complete")drawAll();',
    'else window.addEventListener("load",drawAll);',
    
    '})();',
    '</script>'
  )
}