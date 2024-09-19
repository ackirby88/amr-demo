/**
 * \file   var_tag.h
 * \ingroup amr_group
 * \author akirby
 *
 * \brief AMR tag name information: 
 *      no_tag, box_tag, point_tag, feature_tag, all_tag.
 * 
 * Created on September 7, 2018, 10:46 AM
 */

#ifndef VAR_TAG_H
#define VAR_TAG_H


/** 
 * Name assignment for different amr tags 
 */
enum tag_names {
    no_tag,         /**< No tagging assigned */
    box_tag,        /**< Tagged from box region refinement callback */
    point_tag,      /**< Tagged from point refinement callback */
    feature_tag,    /**< Tagged from feature refinement callback */
    all_tag         /**< Tagged from all refinement callback */
};


#endif /* VAR_TAG_H */

