block_highlight <- function(ggob,tgt_chunk,tgt_block,block_fill="white",border_box="grey70"){
  
  chunks_ys = cbind(seq(from=0, to = 12000, by = 500),seq(from=500, to = 12500, by = 500))
  blocks_xs = cbind(c(.5,4.5,7.5), c(4.5,7.5,9.5))                
  
  myaes_df = data.frame(ymin = chunks_ys[tgt_chunk , 1], ymax = chunks_ys[tgt_chunk , 2],
                      xmin = blocks_xs[tgt_block , 1], xmax = blocks_xs[tgt_block , 2],b_id = paste0("blk_A",LETTERS[1:length(tgt_chunk)]))
  
    outgg = ggob + geom_rect(data=myaes_df, aes(ymin = ymin, ymax = ymax,xmin = xmin, xmax =xmax, group = b_id),
                             inherit.aes = FALSE,alpha=.25, colour=border_box, fill=block_fill)

  return(outgg)
}